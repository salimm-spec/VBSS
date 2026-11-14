clear
clc
close all
%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%%%%%%%
numFFT = 512;
channel_Taps = [1 0.7 0.5 0.35 0.25 0.15 0.1];
channel_Taps = channel_Taps./norm(channel_Taps);
CPL = 20;
bitsPerSubCarrier = 4;                 % 16-QAM
UserCount = [2 3 5 8];
N_usersL=length(UserCount);
qamMapper = comm.RectangularQAMModulator('ModulationOrder',2^bitsPerSubCarrier,'BitInput',true, 'NormalizationMethod','Average power');
qamDemod = comm.RectangularQAMDemodulator( 'ModulationOrder',2^bitsPerSubCarrier,'BitOutput',true,'NormalizationMethod','Average power');
SchemeNames = {'VBSS','FHSS'};
M = 200;         % Total available carriers
N = max(UserCount);          % 5 - 10 Number of users
N_hops = 200;         % Number of hops
qmin =5;        %5 Minimum subband width
qmax =15;        %15 17 - 10  Maximum subband width % min using TaE optimization
qavg=10;          %10 
%%
snr_range = -5:10;
zz=[1];
xx=[50];
for z=1:length(zz)
for x=1:length(xx)
iter_N=zz(z); %% how many times table is repeated befor changed
MC=xx(x); %% mow many indep runs
%%
SNR_L=numel(snr_range);
N_schemes   = numel(SchemeNames);
resultSize = [N_usersL, N_schemes, SNR_L];
BER_Bob       = nan(resultSize);
BER_C_Bob     = nan(resultSize);
BER_Eve       = nan(resultSize);
BER_Eve_conf  = nan(resultSize);
Det_Eve       = nan(resultSize);
Det_Eve_conf  = nan(resultSize);
Asso_Eve      = nan(resultSize);
Asso_Eve_conf = nan(resultSize);
Rec_Eve       = nan(resultSize);
Rec_Eve_conf  = nan(resultSize);
for USL = 1:N_usersL
    N_users = UserCount(USL);
    for schemeIdx = 1:N_schemes
        schemeName = SchemeNames{schemeIdx};
        for snrloop = 1:SNR_L
            snrdb = snr_range(snrloop);
            fprintf('Users=%d | %s | SNR=%d dB | \n', N_users, schemeName, snrdb);
            ber_Leg_Tmp = zeros(1, MC);
            eveRun      = zeros(1, MC);
            detRun      = zeros(1, MC);
            assocRun    = zeros(1, MC);
            recRun      = zeros(1, MC);
            for run_time=1:MC
                rng(run_time+snrloop+5)
                %rng(1)
                if schemeIdx==1
                    %              carr, Users, hops, ......
                    [q,s]=VBSS_Table3(M,N,N_hops,qmin,qmax,qavg);%%%% generate VBSS table
                else
                    [q,s]=VBSS_Table3(M,N,N_hops,10,10,10); %%% generate FHSS Table
                end
                ra = zeros(2, N_hops*iter_N, N_users);
                ra(1,:,:) = repmat(q(:,1:N_users),iter_N,1); % number of carriers assigned to the user
                ra(2,:,:) = repmat(s(:,1:N_users),iter_N,1); % subband starting offset
                %%% data generation for certain SNR, run time
                totalCarUser = sum(ra(1,:,1));
                TDS=totalCarUser*bitsPerSubCarrier*iter_N; % total data size
                gen_data=cell(1,N_users);
                for u = 1:N_users
                    gen_data{u} = randi([0 1],1,TDS);
                end            
                
                bitPtrTx = ones(1,N_users); % current data seg index for transmission
                bitPtrRx = ones(1,N_users); % current data seg index for reception
                rxBits   = zeros(N_users, TDS);
                EveRxAll = cell(1, N_hops * iter_N);
                for hop=1:N_hops*iter_N
                    tableHop = hop;%%% Equivalent hop number
                    txSig=zeros(numFFT+CPL,1);
                    for u=1:N_users
                        nCar = ra(1,tableHop,u); % carriers per user for this hop
                        offset = ra(2,tableHop,u); % subband position for current user in current hop
                        nBits = nCar*bitsPerSubCarrier;
                        idx1 = bitPtrTx(u);
                        idx2 = idx1+nBits-1;
                        seg_data = gen_data{u}(idx1:idx2);
                        bitPtrTx(u) = idx2+1; % increasing index for next hop
                        symbolsIn = step(qamMapper,seg_data');
                        symbolsInOFDM = [zeros(offset,1); symbolsIn; zeros(numFFT-offset-nCar,1)];
                        ifftOut = ifft(symbolsInOFDM,numFFT);
                        txSig = txSig + [ifftOut(numFFT-CPL+1:end);ifftOut]; % CP + Symbols
                    end
                    
                    %%% Legitimate
                    for u = 1:N_users
                        [rxSig,channel_values] = Chan_Ray_AWGN(channel_Taps,txSig,snrdb);
                        RXEQ = RX_module(rxSig,channel_values,numFFT,CPL);
                        nCar = ra(1,tableHop,u);
                        offset = ra(2,tableHop,u);
                        nBits = nCar*bitsPerSubCarrier;
                        RXEQSEL = RXEQ(offset+1:offset+nCar); % extracting the user relative symbols
                        rxbits = step(qamDemod,RXEQSEL);
                        rxbits = rxbits';
                        idx1 = bitPtrRx(u);
                        idx2 = idx1+nBits-1;
                        rxBits(u,idx1:idx2) = rxbits;
                        bitPtrRx(u)=idx2+1; % increasing the index for next hop
                    end
                    %%% Eavesdropper
                    [rxSigEve,channel_valuesEve] = Chan_Ray_AWGN(channel_Taps,txSig,snrdb);
                    EveRxAll{hop} = RX_module(rxSigEve,channel_valuesEve,numFFT,CPL);
                end % end of hops loop
                gen_data_c=zeros(N_users,TDS);
                for u = 1:N_users
                    gen_data_c(u,:) = gen_data{u};
                end 
         err=sum(sum(xor(rxBits,gen_data_c),2));
                ber_Leg_Tmp(run_time)=err/(TDS*N_users);
               [eveBER,eveStats,eveModel] = EavesdropObservation(EveRxAll,N_users,qamDemod,ra,gen_data, bitsPerSubCarrier);
                                                  
                eveRun(run_time) = eveBER;
                detRun(run_time) = eveStats.bandDetectionRate;
                assocRun(run_time) = eveStats.associationAccuracy;
                recRun(run_time) = eveStats.messageRecoveryRate;
            end % end of different run times with diff tables
            BER_Bob(USL,schemeIdx,snrloop)=mean(ber_Leg_Tmp);
            BER_C_Bob(USL,schemeIdx,snrloop)=conf(0.95,ber_Leg_Tmp); % confidence interval
            
            BER_Eve(USL,schemeIdx,snrloop)=mean(eveRun);
            Det_Eve(USL,schemeIdx,snrloop)=mean(detRun);
            Asso_Eve(USL,schemeIdx,snrloop)=mean(assocRun);
            Rec_Eve(USL,schemeIdx,snrloop)=mean(recRun);
            Rec_Eve_C(USL,schemeIdx,snrloop)=conf(0.95,recRun);
            
            %%%% confidence interval
            BER_Eve_conf(USL,schemeIdx,snrloop)=conf(0.95,eveRun);
            Det_Eve_conf(USL,schemeIdx,snrloop)=conf(0.95,detRun);
            Asso_Eve_conf(USL,schemeIdx,snrloop)=conf(0.95,assocRun);
            Rec_Eve_conf(USL,schemeIdx,snrloop)=conf(0.95,recRun);
        end % end of SNR Loop
    end % end of scheme Loop
end% end of user Loop
%%
figure(1);
clf;
% Create main plot
axMain = axes;
hold(axMain,'on');
grid(axMain,'on');
box(axMain,'on');
for u = 1:N_usersL
    Col=['b', 'k'];
    for sc = 1:N_schemes
        y  = squeeze(BER_Bob(u,sc,:));
        ci = squeeze(BER_C_Bob(u,sc,:))/2;
        errorbar(axMain,snr_range,y,ci,'-o', 'LineWidth',1,'MarkerSize',8,'Color',Col(sc));
    end
end
set(axMain,'YScale','log');
xlabel(axMain,'SNR (dB)');
ylabel(axMain,'BER');
title(axMain,'BER Performance');
legend(axMain,SchemeNames,'Location','best');
% Zoomed-in inset
axInset = axes('Position',[0.55 0.55 0.30 0.30]);
hold(axInset,'on');
grid(axInset,'on');
box(axInset,'on');
% Select the user to zoom
uZoom = 1;
Col=['b', 'k'];
for sc = 1:N_schemes
    y  = squeeze(BER_Bob(uZoom,sc,:));
    ci = squeeze(BER_C_Bob(uZoom,sc,:))/2;
    errorbar(axInset,snr_range,y,ci,'-o', 'LineWidth',1,'MarkerSize',8,'Color',Col(sc));
end
set(axInset,'YScale','log');
% Adjust zoom limits
xlim(axInset,[5.5 6.5]);
ylim(axInset,[1e-3 1e-1]);
set(axInset,'FontSize',8);
title(axInset,sprintf('Zoom: %d Users',UserCount(uZoom)));

%% ================= FIGURE 2: EAVESDROPPER BER =================
%figure('Name','Eve BER','Color','w');
figure(2)
hold on;
grid on;
box on;
for u = 1:N_usersL
    Col=['b', 'k'];
    for sc = 1:N_schemes
        y = squeeze(BER_Eve(u,sc,:));
        ci = squeeze(BER_Eve_conf(u,sc,:))/2;
        errorbar(snr_range,y,ci,'-s', 'LineWidth',1,'MarkerSize',6,'Color',Col(sc));
    end      
end
xlabel('SNR (dB)');
ylabel('Eve BER');
legend(SchemeNames,'Location','best');
ylim([0 1]);
title('Eavesdropper BER: VBSS vs FHSS');
%%%% zoomed 
axInset = axes('Position',[0.55 0.55 0.30 0.30]);
hold(axInset,'on');
grid(axInset,'on');
box(axInset,'on');
uZoom = 1;
Col=['b', 'k'];
for sc = 1:N_schemes
    y  = squeeze(BER_Eve(uZoom,sc,:));
    ci = squeeze(BER_Eve_conf(uZoom,sc,:))/2;
    errorbar(axInset,snr_range,y,ci,'-s', 'LineWidth',1,'MarkerSize',8,'Color',Col(sc));
end
set(axInset,'YScale','log');
xlim(axInset,[5.5 6.5]);
ylim(axInset,[.4 .5]);
set(axInset,'FontSize',8);
title(axInset,sprintf('Zoom: %d Users',UserCount(uZoom)));

%% ================= FIGURE 3: EVE ATTACK PERFORMANCE =================
%figure('Name','Eve Attack Performance','Color','w');

for u = 1:N_usersL
figure   
    hold on;
    grid on;
    box on;
Col=['b', 'k'];
    for sc = 1:N_schemes
        plot(snr_range,squeeze(Det_Eve(u,sc,:)),'-o', 'LineWidth',1,'Color',Col(sc));
        plot(snr_range,squeeze(Asso_Eve(u,sc,:)),'--s','LineWidth',1,'Color',Col(sc));
        plot(snr_range,squeeze(Rec_Eve(u,sc,:)),':^','LineWidth',1,'Color',Col(sc));
                   end
    xlabel('SNR (dB)');
    ylabel('Rate');
    title(sprintf('Eavesdropper Attack Performance for %d Users',UserCount(u)));
   legend('VBSS Detection','VBSS Association','VBSS Recovery', ...
       'FHSS Detection','FHSS Association','FHSS Recovery');
    ylim([-0.4 1]);
end
%% ================= FIGURE 4: EAVESDROPPER Rec ================= Rec_Eve
figure()
hold on;
grid on;
box on;
for u = 1:N_usersL
    Col=['b', 'k'];
    for sc = 1:N_schemes
        y = squeeze(Rec_Eve(u,sc,:));
        ci = squeeze(Rec_Eve_C(u,sc,:))/2;
        errorbar(snr_range,y,ci,'-s', 'LineWidth',1,'MarkerSize',8,'Color',Col(sc));
    end      
end
xlabel('SNR (dB)');
ylabel('Recovery Rate');
legend(SchemeNames,'Location','best');
ylim([0 1]);
title('Eavesdropper Recovery Rate: VBSS vs FHSS');
%%%% zoomed 
axInset = axes('Position',[0.55 0.55 0.30 0.30]);
hold(axInset,'on');
grid(axInset,'on');
box(axInset,'on');
uZoom = 2;
Col=['b', 'k'];
for sc = 1:N_schemes
    y  = squeeze(Rec_Eve(uZoom,sc,:));
    ci = squeeze(Rec_Eve_C(uZoom,sc,:))/2;
    errorbar(axInset,snr_range,y,ci,'-s', 'LineWidth',1,'MarkerSize',8,'Color',Col(sc));
end
set(axInset,'YScale','log');
xlim(axInset,[5.5 6.5]);
ylim(axInset,[.1 .2]);
set(axInset,'FontSize',8);
title(axInset,sprintf('Zoom: %d Users',UserCount(uZoom)));

end% x
end %z


