function [rxSig,channel_values]=Chan_Ray_AWGN(channel_Tabs,txSig,snrdb)

 channel_values = channel_Tabs.*(randn(1,7)+1i*randn(1,7))/sqrt(2);
chanOut= conv(channel_values, txSig);
rxSig = awgn(chanOut, snrdb, 'measured'); %% signal in time domain with AWGN

end