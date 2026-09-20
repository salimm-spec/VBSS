%%
clear 
clc
close all
%%
carr=16:48;
Users=[1 2 4];

for nLoop=1:length(carr);
    n=carr(nLoop);
    for NLoop=1:length(Users);
        N=Users(NLoop);
        f=factorial(N);
        combi=nchoosek(n+N,2*N);
        possibilites(nLoop,NLoop)=f*combi;        
        
    end
end
%%
figure
mesh(Users,carr,possibilites)
xlabel('N')
ylabel('n')
grid on