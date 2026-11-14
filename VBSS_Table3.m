function [q,s]=VBSS_Table3(M,N,T,qmin,qmax,qavg) 
%% Initialization
Ctarget=qavg*T;
q = zeros(T,N);       % q(t,i): subband width
s = zeros(T,N);       % s(t,i): starting carrier position
C = zeros(1,N);       % accumulated carriers for each user
R = ones(1,N)*Ctarget;% remainnig carriers for each user
Av_eff=zeros(1,N);
for t=1:T % hops Loop
    userOrder = randperm(N);
    usedCarriers = 0;
    R=Ctarget-C; % remaining =target - been assigned
        % assign width
    for n=1:N
       x=userOrder(n); % current user id 
       if R(x)<=(qmin*(T-t)+qmax)   % 4 
           q(t,x)=qmin;  
       elseif R(x)>=(qmax*(T-t)+qmin)
           q(t,x)=qmax;
       elseif Av_eff(x)>0
           q(t,x)=randi([qmin qavg]);
       elseif Av_eff(x)<0
           q(t,x)=randi([qavg qmax]);
      elseif Av_eff(x)==0
           q(t,x)=randi([qmin qmax]);
      end
        usedCarriers = usedCarriers + q(t,x); % Update number of carriers used in this hop       
        C(x) = C(x) + q(t,x);% Update accumulated number of carriers
    end
   Av_eff=((C/t)-qavg);
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
        maxgap=max(maxgap,0);
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


end