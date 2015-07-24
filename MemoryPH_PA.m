function [PA_Out ScalingGain] = MemoryPH_PA(PAin,PH_f1,PH_f3,PH_f5)

Epsi_1 = PAin;
PH_Branch_1 = filter(PH_f1,1,Epsi_1);
Epsi_3 = PAin.*abs(PAin).^2;
PH_Branch_3 = filter(PH_f3,1,Epsi_3);
Epsi_5 = PAin.*abs(PAin).^4;
PH_Branch_5 = filter(PH_f5,1,Epsi_5);

PA_Out = PH_Branch_1 + PH_Branch_3 + PH_Branch_5;

% Adding noise to the PA model
if 1
    SNR = 70; % dB
    SNR_ratio = 10^(SNR/10);
    Signal_Pwr = mean(abs(PA_Out).^2);
    Noise_pwr = Signal_Pwr/SNR_ratio;
    Noise = sqrt(Noise_pwr)*randn(size(PA_Out));
    PA_Out = PA_Out + Noise;
end
