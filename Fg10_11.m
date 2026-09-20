clear 
clc
close all

%% 
for n=10:512
    FHP(n-9)=1/n;
    VB_SU(n-9)=2/(n*(n+1));
end 

figure
hold on 
semilogy(FHP,'-b')
semilogy(VB_SU,'-r')
legend('FHSS','VBSS')
xlabel('Number of carriers')
ylabel('Probability of detection')
grid on

figure
hold on 
semilogy(1./FHP,'-b')
semilogy(1./VB_SU,'-r')
legend('FHSS','VBSS')
xlabel('Number of carriers')
ylabel('Total posibilities')
grid on
