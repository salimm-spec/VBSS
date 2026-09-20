
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BHSS 5G waveform technique %%%%%
%%%%% Salim A. Mohammed Ali %%%%%%%%%%%%%
%%%%% Edited by KWH %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V4 Evas is function

clear ;
%close all;
clc;

% initializing
numFFT = 512;        % number of FFT points
channel_Tabs = [1 0.7 0.5 0.35 0.25 0.15 0.1];
channel_Tabs=channel_Tabs./norm(channel_Tabs);
CPL=20; % must be > length (channel_Tabs)
SubCarriers=200;
bitsPerSubCarrier=4; % 16 QAM
MC=25; % monte Carlo
qamMapper = comm.RectangularQAMModulator('ModulationOrder',2^bitsPerSubCarrier, 'BitInput', true, 'NormalizationMethod', 'Average power');
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder',2^bitsPerSubCarrier, 'BitOutput', true, 'NormalizationMethod', 'Average power');
% band distribution
NU=10; % for file name
UserCOunt=[2 3 5]; % for simulation
filename=['Lookup_' num2str(NU) '_users.xlsx'];
bandsize = xlsread(filename,'1');
bandPosition = xlsread(filename,'2');

Avg_car_per_hop=10;
snr_range=-5:15; %% SNR Loop
for USL=1:length(UserCOunt)
    
    N_users=UserCOunt(USL) % for simulation
    UsersID=1:N_users;    
    for Loop_rx=1:N_users
        ra_BHSS(1,:,Loop_rx)=bandsize(:,Loop_rx);
        ra_BHSS(2,:,Loop_rx)=bandPosition(:,Loop_rx);
    end    
    seg_no=length(ra_BHSS(1,:,1)); % no of segments    
    total_car=seg_no*Avg_car_per_hop;   
    berb=zeros(N_users,length(snr_range)); % arra for BER calculations
    ber_t=zeros(N_users,length(snr_range)); % total BER of the system
    ber_Ev=zeros(N_users,length(snr_range)); % total BER of the system
    for run_time=1:MC;
        disp(run_time)
        
        %% main Loop
        for snrloop=1:length(snr_range)
            snrdb=snr_range(snrloop);
            disp(snrdb)
            total_data=randi([0 1], N_users, bitsPerSubCarrier*total_car); %% generating total data for three users
            seg_start=ones(N_users,1);
            b_VBSS=zeros(2,total_car*bitsPerSubCarrier,N_users); % bit comparsion array for BHSS
            %%% transmitted bits and received for all segments per one SNR value
            tx_VBSS=[];rx_VBSS=[];EV_VBSS=[];
            for seg_id=1:seg_no
                used_sub_car=sum(ra_BHSS(1,seg_id,:));
                %% %%%% Transmitter
                %        disp('transmission part')
                %%% Proposed BHSS TX
                txSig_BHSS=0;
                tx_d_BHSS=[];
                for Loop_rx=1:N_users
                    seg_data_BHSS=total_data(Loop_rx,seg_start(Loop_rx):ra_BHSS(1,seg_id,Loop_rx)*bitsPerSubCarrier+seg_start(Loop_rx)-1); % selecting tx data from total bits
                    seg_start(Loop_rx)=ra_BHSS(1,seg_id,Loop_rx)*bitsPerSubCarrier+ seg_start(Loop_rx);
                    %                    bsz/ seg no/ user
                    symbolsIn = step(qamMapper,seg_data_BHSS'); % convrting to modulated symbols
                    Start_b=sum(ra_BHSS(1,1:seg_id-1,Loop_rx))*bitsPerSubCarrier+1;
                    End_b=sum(ra_BHSS(1,1:seg_id,Loop_rx))*bitsPerSubCarrier;
                    b_VBSS(1,Start_b:End_b,Loop_rx)=seg_data_BHSS; % 1 refer to tx data % refer to Ln 118
                    tx_d_BHSS=[tx_d_BHSS seg_data_BHSS];
                    offset = ra_BHSS(2,seg_id,Loop_rx);
                    symbolsInOFDM = [zeros(offset,1); symbolsIn;zeros(numFFT-offset-ra_BHSS(1,seg_id,Loop_rx), 1)];
                    ifftOut = ifft(symbolsInOFDM,numFFT);
                    txSig_BHSS=txSig_BHSS+[ifftOut(numFFT-CPL+1:end);ifftOut]; %% adding cyclic prefix
                end
                %%% End of transmission   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                %% %%%%%%       Channel and Receiver
                %%%%%%%%%% BHSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bit_per_seg=0;
                rxBits_seg=[];%zeros(1,used_sub_car*bitsPerSubCarrier); % refer to Ln 71
                for Loop_rx=1:N_users+1
                    if Loop_rx<N_users+1
                        %%% Users
                        [rxSig,channel_values]=Chan_Ray_AWGN(channel_Tabs,txSig_BHSS,snrdb);
                        RXEQ=RX_module(rxSig,channel_values,numFFT,CPL);
                        offset = ra_BHSS(2,seg_id,Loop_rx);
                        RXEQSEL=RXEQ(offset+(1:ra_BHSS(1,seg_id,Loop_rx))); %% select in band data
                        rxBits = step(qamDemod,RXEQSEL);
                        Start_b=sum(ra_BHSS(1,1:seg_id-1,Loop_rx))*bitsPerSubCarrier+1;
                        End_b=sum(ra_BHSS(1,1:seg_id,Loop_rx))*bitsPerSubCarrier;
                        b_VBSS(2,Start_b:End_b,Loop_rx)=rxBits';
                        rxBits_seg=[rxBits_seg rxBits']; % received bits y all users for current segment
                    else
                        %%% Eavesdroper %%%%%%%%%%%%%%%
                        [rxSig,channel_values]=Chan_Ray_AWGN(channel_Tabs,txSig_BHSS,snrdb); %% Eve Channel
                        RXEQ_EVE_BHSS=RX_module(rxSig,channel_values,numFFT,CPL); %% FFT +CP removal
                        [Ev_thr_BHSS,Ev_topK_BHSS]=Eavesdrop(RXEQ_EVE_BHSS,qamDemod,used_sub_car);
                    end
                end
                tx_VBSS=[tx_VBSS tx_d_BHSS];
                rx_VBSS=[rx_VBSS rxBits_seg];
                EV_VBSS=[EV_VBSS Ev_topK_BHSS'];
                %  disp(seg_id)
            end  %%% Segment Loop
            err_VBSS_S_Lg(snrloop,run_time)=sum(xor(tx_VBSS,rx_VBSS));
            err_VBSS_S_EV(snrloop,run_time)=sum(xor(tx_VBSS,EV_VBSS));
            
            for u_id=1:N_users
                err_VBSS_U_Lg(snrloop,run_time,u_id)= sum(xor(b_VBSS(1,:,u_id),b_VBSS(2,:,u_id)));
            end
        end     %% SNR Loop
        
        
    end % monte carlo loop (run_time)
    %%
    dume1=run_time*N_users*total_car*bitsPerSubCarrier;
    dume2=run_time*total_car*bitsPerSubCarrier;
    %% BER for whole system per method
    ber_VBSS_S_Lg=sum(err_VBSS_S_Lg,2)/dume1;
    ber_VBSS_S_EV=sum(err_VBSS_S_EV,2)/dume1;
    % BER for indiviual users
    ber_VBSS_U_Lg=permute(sum(err_VBSS_U_Lg,2)/dume2, [1 3 2]);
    %%
    xlswrite('VBSS_SY_Lg.xlsx',ber_VBSS_S_Lg,N_users)
    
    %%
    figure(11)
    hold on
    Legnd=[];
    for Loop_rx=1:N_users
        semilogy(snr_range,ber_VBSS_U_Lg(:,Loop_rx),'LineWidth',2)
        Legnd=[Legnd , ['user ' num2str(Loop_rx)]];
    end
    
    legend('user1','user2','user3')
    xlabel('snr (dB)')
    ylabel('BER')
    %%
    figure(12)
    hold on
    semilogy(snr_range,ber_VBSS_S_Lg,'m-.','LineWidth',2)
    semilogy(snr_range,ber_VBSS_S_EV,'m-.','LineWidth',2)
    
    legend('Legitimate','Eavesdropper')
    xlabel('SNR (dB)')
    ylabel('BER')
end

