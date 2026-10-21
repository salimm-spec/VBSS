
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
Lee=[];
filename = ['Lookup_' num2str(NU) '_users.xlsx'];
bandsize = xlsread(filename,'1');
bandPosition = xlsread(filename,'2');
FHSS_tables = cell(1,length(UserCOunt));
for USL = 1:length(UserCOunt)
    N_users = UserCOunt(USL);
    FHSS_tables{USL} = makeFHSSAllocation( numFFT,N_Hops,N_users,FHSS_bandwidth,FHSS_seed+N_users);
end
Result = struct();% Results are stored separately for the two schemes.
%% ===================
rng(1)
MC = 3;                                % use 25 for final results
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
                fprintf('Monte-Carlo run %d/%d  SNR = %d dB\n',run_time,MC,snrdb);
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
                    txAll = [txAll txHopBits]; %#ok<AGROW>
                    rxAll = [rxAll rxHopBits]; %#ok<AGROW>
                end               
                err_Lg(snrloop,run_time) = sum(xor(txAll,rxAll)); %  LEGITIMATE BER 
                
                
            end
        end
        %========================= AVERAGES =============================
        ber_Lg = mean(err_Lg,2)/totalBits;
 
        Result.(schemeName)(USL).N_users = N_users;
        Result.(schemeName)(USL).snr = snr_range;
        Result.(schemeName)(USL).legitimateBER = ber_Lg;

        %% =========================== PLOTS =============================
        figure(1000); 
        semilogy(snr_range,ber_Lg,'LineWidth',2); hold on;
        grid on;
        xlabel('SNR (dB)');
        ylabel('BER');
       Lee=[Lee; num2str(N_users) ' Legitimate Users for ' schemeName];
               
    end    
end
 legend(Lee,'Location','southwest');
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

    % Random permutation of possible starting positions. Select positions
    % greedily while maintaining non-overlapping fixed-width bands.
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

    % Randomize user-to-band assignment while preserving the same set of
    % non-overlapping bands.
    starts = starts(randperm(N_users));

    for u = 1:N_users
        ra_FHSS(1,hop,u) = q;
        ra_FHSS(2,hop,u) = starts(u);
    end
end
end
