function [eveBER,stats,model] = EavesdropObservation( ...
    EveRxAll,N_users,qamDemod,ra,total_data,bitsPerSubCarrier)



Hops = length(EveRxAll);
% 1. BLIND ENERGY DETECTION
obs = struct('hop',{},'start',{},'stop',{},'width',{},'energy',{},'energyWidth',{},'bits',{});
nObs = 0;
for hop = 1:Hops
    Y = EveRxAll{hop};
    P = abs(Y(:)).^2;
    Psm = movmean(P,3);
    % Robust noise-floor estimate. No allocation-table information.
    noiseMedian = median(Psm);
    noiseMAD = median(abs(Psm-noiseMedian));
    noiseSigma = 1.4826*noiseMAD + eps;
    threshold = noiseMedian + 4*noiseSigma;
    occupied = Psm > threshold;
    % Fill small holes inside an occupied contiguous band.
    occupied = conv(double(occupied),ones(3,1),'same') >= 2;
    edge = diff([false;occupied;false]);
    starts = find(edge==1);
    stops  = find(edge==-1)-1;

    keep = false(length(starts),1);
    for k = 1:length(starts)
        idx = starts(k):stops(k);
        E = mean(Psm(idx));
        W = length(idx);
        keep(k) = (E > noiseMedian + 2*noiseSigma) && (W >= 2);
    end
    starts = starts(keep);
    stops = stops(keep);
    for k = 1:length(starts)
        idx = starts(k):stops(k);
        W = length(idx);
        E = mean(Psm(idx));
        nObs = nObs + 1;
        obs(nObs).hop = hop;
        obs(nObs).start = starts(k);
        obs(nObs).stop = stops(k);
        obs(nObs).width = W;
        obs(nObs).energy = E;
        obs(nObs).energyWidth = E/max(W,1);
        % Demodulation is performed without assigning a user label.
        b = step(qamDemod,Y(idx));
        obs(nObs).bits = b(:)';
    end
end

model = struct;
model.valid = false;
model.Hops = Hops;
model.N_users = N_users;
model.observations = obs;
model.numObservations = nObs;

stats = struct('bandDetectionRate',0, ...
    'associationAccuracy',0, ...
    'conditionalAssociationAccuracy',0, ...
    'messageRecoveryRate',0, ...
    'attemptedBER',NaN, ...
    'detectedBandCount',0, ...
    'trueBandCount',Hops*N_users, ...
    'correctAssociationCount',0, ...
    'attemptedBits',0, ...
    'recoveredCorrectBits',0, ...
    'bitErrors',0);

totalTrueBits = 0;
for u = 1:N_users
    totalTrueBits = totalTrueBits + ...
        sum(ra(1,:,u))*bitsPerSubCarrier;
end

if nObs < N_users
    model.failureReason = 'Too few blind detections for K-means.';
    eveBER = NaN;
    return;
end
% 2. OBSERVABLE FEATURES
F = zeros(nObs,3);
for n = 1:nObs
    F(n,1) = log(1 + obs(n).width);
    F(n,2) = log(1 + obs(n).energy);
    F(n,3) = log(1 + obs(n).energyWidth);
end

mu = mean(F,1);
sigma = std(F,0,1);
sigma(sigma < 1e-12) = 1;
Z = (F-mu)./sigma;

% 3. K-MEANS OVER ALL OBSERVATIONS
rng(20261004,'twister');
opts = statset('MaxIter',1000,'Display','off');

[clusterID,centers] = kmeans(Z,N_users, ...
    'Distance','sqeuclidean', ...
    'Replicates',50, ...
    'Options',opts);

model.features = F;
model.standardizedFeatures = Z;
model.clusterID = clusterID;
model.clusterCenters = centers;
model.featureMean = mu;
model.featureStd = sigma;

model.anonymousGroups = cell(1,N_users);
for c = 1:N_users
    model.anonymousGroups{c} = find(clusterID==c);
end

% 4. POST-ATTACK MATCHING FOR EVALUATION ONLY

matches = [];  % [obsIndex trueUser IoU]

for hop = 1:Hops

    detIDs = find(arrayfun(@(x) x.hop==hop,obs));
    if isempty(detIDs)
        continue;
    end

    candidates = [];

    for d = 1:length(detIDs)
        n = detIDs(d);
        ds = obs(n).start;
        de = obs(n).stop;

        for u = 1:N_users
            ts = ra(2,hop,u)+1;
            te = ts + ra(1,hop,u)-1;

            inter = max(0,min(de,te)-max(ds,ts)+1);
            unionLen = (de-ds+1) + (te-ts+1) - inter;
            iou = inter/max(unionLen,1);

            if iou >= 0.50
                candidates = [candidates; n u iou]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(candidates)
        % Strongest spatial matches first.
        [~,ord] = sort(candidates(:,3),'descend');
        candidates = candidates(ord,:);

        usedDet = false(length(detIDs),1);
        usedUser = false(1,N_users);

        for k = 1:size(candidates,1)
            n = candidates(k,1);
            u = candidates(k,2);
            detLocal = find(detIDs==n,1);

            if ~usedDet(detLocal) && ~usedUser(u)
                matches = [matches; candidates(k,:)]; %#ok<AGROW>
                usedDet(detLocal) = true;
                usedUser(u) = true;
            end
        end
    end
end

model.matches = matches;
model.detectedBandCount = size(matches,1);

% 5. DETECTION METRIC
trueBandCount = Hops*N_users;
detectionRate = size(matches,1)/max(trueBandCount,1);

% 6. LEARN ANONYMOUS CLUSTER -> USER ASSOCIATION

if size(matches,1) == 0
    eveBER = NaN;
    model.failureReason = 'No valid detected bands.';
    return;
end

associationScore = zeros(N_users,N_users);

for k = 1:size(matches,1)
    n = matches(k,1);
    trueUser = matches(k,2);
    c = clusterID(n);
    associationScore(c,trueUser) = associationScore(c,trueUser)+1;
end

P = perms(1:N_users);
bestScore = -inf;
bestMap = 1:N_users;

for rowIdx = 1:size(P,1)
    s = 0;
    
    for c = 1:N_users
        s = s + associationScore(c,P(rowIdx, c));
    end
    if s > bestScore
        bestScore = s;
        bestMap = P(rowIdx, :);
    end
end

model.associationScore = associationScore;
model.clusterToUser = bestMap;

% 7. ASSOCIATION METRICS
correctAssociationCount = 0;

for k = 1:size(matches,1)
    n = matches(k,1);
    trueUser = matches(k,2);
    estimatedUser = bestMap(clusterID(n));

    if estimatedUser == trueUser
        correctAssociationCount = correctAssociationCount + 1;
    end
end
associationAccuracy = correctAssociationCount/max(trueBandCount,1);
conditionalAssociationAccuracy = correctAssociationCount/max(size(matches,1),1);
% 8. LEARNED ASSOCIATION RECORD
association = zeros(size(matches,1),8);

for k = 1:size(matches,1)
    n = matches(k,1);
    trueUser = matches(k,2);
    iou = matches(k,3);
    estimatedUser = bestMap(clusterID(n));

    association(k,:) = [obs(n).hop,obs(n).start,obs(n).stop, clusterID(n),estimatedUser,trueUser,iou,estimatedUser==trueUser];
end
model.association = association;
% 9. 16-QAM MESSAGE RECOVERY
attemptedBits = 0;
recoveredCorrectBits = 0;
bitErrors = 0;
recoveryRecord = [];

for k = 1:size(matches,1)

    n = matches(k,1);
    trueUser = matches(k,2);
    estimatedUser = bestMap(clusterID(n));

    if estimatedUser ~= trueUser
        continue;
    end

    hop = obs(n).hop;
    nCar = ra(1,hop,trueUser);
    nBits = nCar*bitsPerSubCarrier;

    bitStart = 1 + sum(ra(1,1:hop-1,trueUser))*bitsPerSubCarrier;
    bitStop = bitStart+nBits-1;
    trueBits = total_data{trueUser}(bitStart:bitStop);

    recovered = obs(n).bits;
    nCompare = min(length(recovered),nBits);

    if nCompare > 0
        e = sum(xor(trueBits(1:nCompare),recovered(1:nCompare)));
        bitErrors = bitErrors + e;
        attemptedBits = attemptedBits + nCompare;
        recoveredCorrectBits = recoveredCorrectBits + (nCompare-e);
    end

    recoveryRecord = [recoveryRecord;hop trueUser n nCar nBits nCompare]; %#ok<AGROW>
end

messageRecoveryRate = recoveredCorrectBits/max(totalTrueBits,1);

if attemptedBits > 0
    eveBER = bitErrors/attemptedBits;
else
    eveBER = NaN;
end

stats.bandDetectionRate = detectionRate;
stats.associationAccuracy = associationAccuracy;
stats.conditionalAssociationAccuracy = conditionalAssociationAccuracy;
stats.messageRecoveryRate = messageRecoveryRate;
stats.attemptedBER = eveBER;
stats.detectedBandCount = size(matches,1);
stats.trueBandCount = trueBandCount;
stats.correctAssociationCount = correctAssociationCount;
stats.attemptedBits = attemptedBits;
stats.recoveredCorrectBits = recoveredCorrectBits;
stats.bitErrors = bitErrors;

model.trueBandCount = trueBandCount;
model.detectedBandCount = size(matches,1);
model.correctAssociationCount = correctAssociationCount;
model.associationAccuracy = associationAccuracy;
model.conditionalAssociationAccuracy = conditionalAssociationAccuracy;
model.totalTrueBits = totalTrueBits;
model.attemptedBits = attemptedBits;
model.recoveredCorrectBits = recoveredCorrectBits;
model.bitErrors = bitErrors;
model.recoveryRecord = recoveryRecord;
model.bestPermutationScore = bestScore;
model.valid = true;

clusterSummary = zeros(N_users,5);
for c = 1:N_users
    ids = find(clusterID==c);
    if isempty(ids)
        continue;
    end

    w = zeros(length(ids),1);
    e = zeros(length(ids),1);
    ew = zeros(length(ids),1);

    for j = 1:length(ids)
        w(j) = obs(ids(j)).width;
        e(j) = obs(ids(j)).energy;
        ew(j) = obs(ids(j)).energyWidth;
    end

    clusterSummary(c,:) = [mean(w),mean(e),mean(ew),length(ids),bestMap(c)];
end
model.clusterSummary = clusterSummary;

end
