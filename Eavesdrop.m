function [Ev_thr,Ev_topK]=Eavesdrop(RXEQ_EVE,qamDemod,Used_Car)

%%% energy threshold based
RXEQ_energy=abs(RXEQ_EVE).^2;
power_thresh = mean(RXEQ_energy);   % or k * noise_variance  % or 1
RXEQ_pos_thr=find(sign(RXEQ_energy-power_thresh)+1);
RXEQSEL_thr=RXEQ_EVE(RXEQ_pos_thr); %% only those pass threshold
Ev_thr = step(qamDemod,RXEQSEL_thr);




%%% top k energy based detection
[~,RXEQ_pos]=sort(RXEQ_energy,'descend');
RXEQ_pos_ord=sort(RXEQ_pos(1:Used_Car));
RXEQSEL=RXEQ_EVE(RXEQ_pos_ord);
Ev_topK = step(qamDemod,RXEQSEL);
end