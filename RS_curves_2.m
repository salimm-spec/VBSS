%%
clc
clear
close all

n=[128 200 265 512];
N=[3 5 10 20 50 100];
SNR=-5:5:35;

Rs1=zeros(length(n),length(N),length(SNR));
Rs2=zeros(length(n),length(N),length(SNR));

for car_loop=1:length(n)
    nv=n(car_loop); %value of carriers
    for user_loop=1:length(N)
        Nv=N(user_loop); % value of users
        P1=1/Nv;
        P2=1/(factorial(Nv));
        for SNR_loop=1:length(SNR);
            SNRv=10^(SNR(SNR_loop)/10);           
            
            nu=1+SNRv;
            de_1=1+P1*SNRv;
            de_2=1+P2*SNRv;
            
            Rs1(car_loop,user_loop,SNR_loop)=((nv+Nv)*log2(nu/de_1))/2;
            Rs2(car_loop,user_loop,SNR_loop)=((nv+Nv)*log2(nu/de_2))/2;
        end
    end
end
RS(:,:,:,1)=Rs1;
RS(:,:,:,2)=Rs2;
pc=0;
for i=1:2
    Rs=RS(:,:,:,i);

    %% fig 1
figure ()
R=permute(Rs,[3 2 1]);
plot(SNR,R(:,:,2),'LineWidth',2)
xlabel('SNR (dB)','FontSize', 14)
ylabel('R_s','FontSize', 14)
legend('N=3', 'N=5', 'N=10', 'N=20', 'N=50', 'N=100')
title(['Secrecy rate for n=200 and Rs' num2str(i)])
grid on
box on 
%% fig 2
% SNR=20 dB Rs(:,:,6)
figure ()
R=Rs;
plot(N,R(:,:,2),'LineWidth',2)
xlabel('number of users','FontSize', 14)
ylabel('R_s','FontSize', 14)
legend('n=128', 'n=200', 'n=265', 'n=512')
title(['Secrecy rate for SNR = 20 dB and Rs' num2str(i)])
grid on
box on 
%% fig 3
figure ()
R=permute(Rs,[3 2 1]);
mesh(N,SNR,R(:,:,2))
xlabel('number of users','FontSize', 14)
ylabel('SNR (dB)','FontSize', 14)
zlabel('R_s','FontSize', 14)
title(['Secrecy rate for SNR = 20 dB and Rs' num2str(i)])
grid on
box on 
end
