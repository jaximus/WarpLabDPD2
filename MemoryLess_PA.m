function PA_Out = MemoryLess_PA(PAin,Beta_1,Beta_3,Beta_5)

Epsi_1 = PAin;
Epsi_3 = PAin.*abs(PAin).^2;
Epsi_5 = PAin.*abs(PAin).^4;

PA_Out = Beta_1*Epsi_1 + Beta_3*Epsi_3 + Beta_5*Epsi_5;

% Adding noise to the PA model
if 1
    SNR = 70; % dB
    SNR_ratio = 10^(SNR/10);
    Signal_Pwr = mean(abs(PA_Out).^2);
    Noise_pwr = Signal_Pwr/SNR_ratio;
    Noise = sqrt(Noise_pwr)*randn(size(PA_Out));
    PA_Out = PA_Out + Noise;
end


