
clear;
clc;
close all;
numFFT = 512;
channel_Taps = [1 0.7 0.5 0.35 0.25 0.15 0.1];
channel_Taps = channel_Taps./norm(channel_Taps);
CPL = 20;
bitsPerSubCarrier = 4;                 % 16-QAM
UserCOunt = [2 3 5];
N_Hops = 200;
FHSS_bandwidth = 10;
FHSS_seed = 73129; % Reproducible FHSS hopping table.
qamMapper = comm.RectangularQAMModulator('ModulationOrder',2^bitsPerSubCarrier,'BitInput',true, 'NormalizationMethod','Average power');
qamDemod = comm.RectangularQAMDemodulator( 'ModulationOrder',2^bitsPerSubCarrier,'BitOutput',true,'NormalizationMethod','Average power');
NU = 5;
filename = ['Lookup_' num2str(NU) '_users.xlsx'];
bandsize = xlsread(filename,'1');
bandPosition = xlsread(filename,'2');
FHSS_tables = cell(1,length(UserCOunt));
for USL = 1:length(UserCOunt)
    N_users = UserCOunt(USL);
    FHSS_tables{USL} = makeFHSSAllocation(numFFT,N_Hops,N_users,FHSS_bandwidth,FHSS_seed+N_users);
end
Result = struct();% Results are stored separately for the two schemes.
%% ==============================
rng(1)
MC = 100; 
snr_range = -5:15;
%% ============================ SIMULATION ===============================
for USL = 1:length(UserCOunt)    
    N_users = UserCOunt(USL);
    ra_VBSS = zeros(2,N_Hops,N_users);
    for u = 1:N_users
        ra_VBSS(1,:,u) = bandsize(1:N_Hops,u);
        ra_VBSS(2,:,u) = bandPosition(1:N_Hops,u);
    end
    ra_FHSS = FHSS_tables{USL};
    SchemeNames = {'VBSS','FHSS'};
    SchemeTables = {ra_VBSS,ra_FHSS};
    for schemeIdx = 1:2
        schemeName = SchemeNames{schemeIdx};
        ra = SchemeTables{schemeIdx};
        totalCarUser = zeros(1,N_users);
        for u = 1:N_users
            totalCarUser(u) = sum(ra(1,:,u));
        end
       
        totalBits = sum(totalCarUser)*bitsPerSubCarrier;
        
        err_Lg = zeros(length(snr_range),MC);
        err_Eve = nan(length(snr_range),MC);
        
        eve_detection = zeros(length(snr_range),MC);
        eve_association = zeros(length(snr_range),MC);
        eve_conditionalAssociation = zeros(length(snr_range),MC);
        eve_recovery = zeros(length(snr_range),MC);
        eve_attempted_ber = nan(length(snr_range),MC);
        % ======================= MONTE CARLO ============================
        for run_time = 1:MC
            for snrloop = 1:length(snr_range)
                snrdb = snr_range(snrloop);
                fprintf('Monte-Carlo run %d/%d  SNR = %d dB %s %d users\n',run_time,MC,snrdb,schemeName,N_users);
                total_data = cell(1,N_users);
                for u = 1:N_users
                    total_data{u} = randi([0 1],1,totalCarUser(u)*bitsPerSubCarrier);
                end
                bitPtr = ones(1,N_users);
                txAll = []; rxAll = []; EveRxAll = cell(1,N_Hops);
                for hop = 1:N_Hops
                    txSig = 0;
                    txHopBits = [];
                    for u = 1:N_users
                        nCar = ra(1,hop,u);
                        offset = ra(2,hop,u);
                        nBits = nCar*bitsPerSubCarrier;
                        if offset < 0 || offset+nCar > numFFT
                            error('%s allocation exceeds FFT size at hop %d, user %d.',schemeName,hop,u);
                        end
                        idx1 = bitPtr(u);
                        idx2 = bitPtr(u)+nBits-1;
                        seg_data = total_data{u}(idx1:idx2);
                        bitPtr(u) = idx2+1;
                        symbolsIn = step(qamMapper,seg_data');
                        txHopBits = [txHopBits seg_data]; %#ok<AGROW>
                        symbolsInOFDM = [zeros(offset,1); symbolsIn; zeros(numFFT-offset-nCar,1)];
                        ifftOut = ifft(symbolsInOFDM,numFFT);
                        txSig = txSig + [ifftOut(numFFT-CPL+1:end);ifftOut];
                    end
                    %---------------- LEGITIMATE RECEIVERS --------------
                    rxHopBits = [];
                    for u = 1:N_users
                        [rxSig,channel_values] = Chan_Ray_AWGN(channel_Taps,txSig,snrdb);
                        RXEQ = RX_module(rxSig,channel_values,numFFT,CPL);
                        nCar = ra(1,hop,u);
                        offset = ra(2,hop,u);
                        RXEQSEL = RXEQ(offset+(1:nCar));
                        rxBits = step(qamDemod,RXEQSEL);
                        rxBits = rxBits(:)';
                        rxHopBits = [rxHopBits rxBits]; %#ok<AGROW>
                    end
                    
                    % ---------------------- EVE ------------------------
                    [rxSigEve,channel_valuesEve] = Chan_Ray_AWGN(channel_Taps,txSig,snrdb);
                    RXEQ_EVE = RX_module(rxSigEve,channel_valuesEve,numFFT,CPL);
                    EveRxAll{hop} = RXEQ_EVE;
                    txAll = [txAll txHopBits]; %#ok<AGROW>
                    rxAll = [rxAll rxHopBits]; %#ok<AGROW>
                end
                err_Lg(snrloop,run_time) = sum(xor(txAll,rxAll)); %  LEGITIMATE BER
                see=rng;
                [eveBER,eveStats,eveModel] = EavesdropObservation( EveRxAll,N_users,qamDemod,ra,total_data, bitsPerSubCarrier);
                rng(see)
                err_Eve(snrloop,run_time) = eveBER;
                eve_detection(snrloop,run_time) = eveStats.bandDetectionRate;
                eve_association(snrloop,run_time) = eveStats.associationAccuracy;
                eve_conditionalAssociation(snrloop,run_time) = eveStats.conditionalAssociationAccuracy;
                eve_recovery(snrloop,run_time) = eveStats.messageRecoveryRate;
                eve_attempted_ber(snrloop,run_time) = eveStats.attemptedBER;
%                 if run_time == MC
%                     save(sprintf('EveModel_%s_%dUsers_SNR%d.mat', schemeName,N_users,snrdb),'eveModel');
%                 end
            end
        end
        %========================= AVERAGES =============================
        ber_Lg = mean(err_Lg,2)/totalBits;
        ber_Eve = nan(length(snr_range),1);
        for ii = 1:length(snr_range)
            vals = err_Eve(ii,~isnan(err_Eve(ii,:)));
            if ~isempty(vals)
                ber_Eve(ii) = mean(vals);
            end
        end
        meanDetection = mean(eve_detection,2);
        meanAssociation = mean(eve_association,2);
        meanConditionalAssociation = mean(eve_conditionalAssociation,2);
        meanRecovery = mean(eve_recovery,2);
        Result.(schemeName)(USL).N_users = N_users;
        Result.(schemeName)(USL).snr = snr_range;
        Result.(schemeName)(USL).legitimateBER = ber_Lg;
        Result.(schemeName)(USL).eveBER = ber_Eve;
        Result.(schemeName)(USL).detection = meanDetection;
        Result.(schemeName)(USL).association = meanAssociation;
        Result.(schemeName)(USL).conditionalAssociation = meanConditionalAssociation;
        Result.(schemeName)(USL).recovery = meanRecovery;
        
        
    end
    % Direct comparison plots
    figure(300+USL); clf;
    plot(snr_range,Result.VBSS(USL).detection,'LineWidth',2,'Color','b','Marker','+'); hold on;
    plot(snr_range,Result.FHSS(USL).detection,'--','LineWidth',2,'Color','k','Marker','s');
    plot(snr_range,Result.VBSS(USL).association,'LineWidth',2,'Color','b','Marker','+');
    plot(snr_range,Result.FHSS(USL).association,'--','LineWidth',2,'Color','k','Marker','s');
    plot(snr_range,Result.VBSS(USL).recovery,'LineWidth',2,'Color','b','Marker','+');
    plot(snr_range,Result.FHSS(USL).recovery,'--','LineWidth',2,'Color','k','Marker','s');
    grid on;
    ylim([0 1]);
    xlabel('SNR (dB)');
    ylabel('Rate');
    title(sprintf('VBSS vs FHSS - Observation-based Eve (%d users)',N_users));
    legend('VBSS Detection','FHSS Detection','VBSS Association','FHSS Association','VBSS Recovery','FHSS Recovery', 'Location','southeast');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(100+USL);
    semilogy(snr_range,Result.VBSS(USL).legitimateBER,'LineWidth',2,'Color','b','Marker','+'); hold on;
    semilogy(snr_range,Result.VBSS(USL).eveBER,'LineWidth',2,'Color','b','Marker','s');
    semilogy(snr_range,Result.FHSS(USL).legitimateBER,'LineWidth',2,'Color','k','Marker','+');
    semilogy(snr_range,Result.FHSS(USL).eveBER,'LineWidth',2,'Color','k','Marker','s');
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title(sprintf('%d-user : Legitimate vs Observation-based Eve',N_users));
    legend('Legitimate_VBSS','Eve_VBSS','Legitimate_FHSS','Eve_FHSS','Location','southwest');
    
end

save('VBSS_FHSS_Observation_Results.mat','Result');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTION: FIXED-BAND FHSS ALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ra_FHSS = makeFHSSAllocation(numFFT,N_Hops,N_users,q,seed)
if N_users*q > numFFT
    error('FHSS allocation impossible: N_users*q exceeds numFFT.');
end

ra_FHSS = zeros(2,N_Hops,N_users);
rng(seed,'twister');

for hop = 1:N_Hops
    possible = randperm(numFFT-q+1);
    starts = zeros(1,N_users);
    nSelected = 0;
    for ii = 1:length(possible)
        s = possible(ii)-1;              % zero-based start
        if nSelected == 0
            ok = true;
        else
            ok = true;
            for jj = 1:nSelected
                oldS = starts(jj);
                if ~(s+q <= oldS || oldS+q <= s)
                    ok = false;
                    break;
                end
            end
        end
        
        if ok
            nSelected = nSelected + 1;
            starts(nSelected) = s;
            if nSelected == N_users
                break;
            end
        end
    end
    
    if nSelected < N_users
        error('Could not create a non-overlapping FHSS allocation at hop %d.',hop);
    end
     starts = starts(randperm(N_users));
    
    for u = 1:N_users
        ra_FHSS(1,hop,u) = q;
        ra_FHSS(2,hop,u) = starts(u);
    end
end
end
