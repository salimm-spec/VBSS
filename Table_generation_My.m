clc;
clear;
close all;
M = 200;         % Total available carriers
N = 5;          % 5 - 10 Number of users
T = 200;         % Number of hops
qmin = 5;        % Minimum subband width
qmax =15;        % 17 - 10  Maximum subband width % min using TaE optimization
qavg=10;          % 10 - 7
saftygap= 4;% mod(M/N,qmax-qmin)+1 ;      % 4 - 2 is optimized using trial and error
Ctarget=qavg*T;
q = zeros(T,N);       % q(t,i): subband width
s = zeros(T,N);       % s(t,i): starting carrier position
C = zeros(1,N);       % accumulated carriers for each user
R = ones(1,N)*Ctarget;% remainnig carriers for each user

for t=1:T % hops Loop
    userOrder = randperm(N);
    usedCarriers = 0;
    R=Ctarget-C; % remaining =target - been assigned
    % assign width
    for n=1:N
       x=userOrder(n); % current user id 
       min_required=(N-n)*qmin; 
       max_allowed= min(qmax, R(x));
       max_allowed= min(max_allowed, M-usedCarriers-min_required);       
       if R(x)<=qmin*(T-t+saftygap)   % 4 
           q(t,x)=qmin; %randi([qmin max_allowed]); 
       else
            q(t,x)=randi([qmin max_allowed]);
       end
        usedCarriers = usedCarriers + q(t,x); % Update number of carriers used in this hop       
        C(x) = C(x) + q(t,x);% Update accumulated number of carriers
       
    end
end
for n=1:N
  if Ctarget-C(n)> 0 
     q(T,n)=q(T,n)+(Ctarget-C(n)) ;  
  end 
end
C=sum(q,1);
for t=1:T
   currentPosition = randi([0 5]);
   userOrder = randperm(N);
   maxgap=M-sum(q(t,:),2)-currentPosition;
    for n = 1:N 
        x=userOrder(n);
        currentgap=randi([0 maxgap]);
        maxgap=maxgap-currentgap;
        if q(t,x) > 0            
            s(t,x) = currentPosition;            
            currentPosition = currentPosition + q(t,x)+ currentgap;           
        else
            s(t,x) = NaN;            
        end        
    end  
    
end

filename=['Lookup_' num2str(N) '_users.xlsx'];
xlswrite(filename,q,'1')
xlswrite(filename,s,'2')
