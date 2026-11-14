function C_val=conf(ar,ber_Leg_Tmp)
L=length(ber_Leg_Tmp);
crit=tinv(ar,L-1);
Stdv=std(ber_Leg_Tmp);
C_val=2*crit*Stdv/L;
end