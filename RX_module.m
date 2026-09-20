
function RXEQ=RX_module(rxSig,channel_values,numFFT,CPL)
Rx = rxSig(CPL+1:CPL+numFFT);
RxSymbols = fft(Rx,numFFT); %% signal in frequency domain , should equal to OFDMS_CK
chanF=(fft(channel_values,numFFT)).'; %% channel equvalent
RXEQ=RxSymbols ./ chanF;
end