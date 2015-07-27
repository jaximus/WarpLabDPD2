function LoopDelay = DPD_LoopDelayEst(PA_InputSignal, IM3GeneratedSignal, MemoryLessPA, SystemFs,Signal_Bandwidth, ...
                             IM3_Freq,Beta_1,Beta_3,Beta_5,PH_f1,PH_f3,PH_f5,AdditionalDelay)

% Third order IMD extraction filter, all frequency values are in MHz.
Fs    = round(SystemFs/1e6);  % Sampling Frequency
N     = 200;  % Order
Fpass = (3*Signal_Bandwidth)/2; % Passband Frequency
Fstop = (5*Signal_Bandwidth)/2; % Stopband Frequency
Wpass = 1;    % Passband Weight
Wstop = 5;    % Stopband Weight
dens  = 20;   % Density Factor
% Calculate the coefficients using the FIRPM function.
IM3Filter  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});

% Reset the delay line
DelayLineReset = 1;
DelayLine(0,0,DelayLineReset);
DelayLineReset = 0;

NumSamples = 10000;
ErrorSignal = zeros(1,NumSamples);

for Sample = 1:NumSamples
    
    PA_InSample = PA_InputSignal(Sample);
    
    if MemoryLessPA
        PA_OutSample = MemoryLess_SampleBased_PA(PA_InSample,Beta_1,Beta_3,Beta_5);
    else
        PA_OutSample = MemoryPH_SampleBased_PA(PA_InSample,PH_f1,PH_f3,PH_f5);
    end
    
    PA_OutSampleShifted = PA_OutSample*exp(-2*pi*1i*Sample*IM3_Freq*1e6/SystemFs);
    
    if Sample == 1
        [IM3FilteredSample Current_State] = filter(IM3Filter,1,PA_OutSampleShifted);
    else
        [IM3FilteredSample Current_State] = ... 
            filter(IM3Filter,1,PA_OutSampleShifted,Previous_State);
    end
    Previous_State = Current_State;
    
    if AdditionalDelay ~= 0
        IM3FilteredSampleDelayed = DelayLine(IM3FilteredSample,AdditionalDelay,DelayLineReset);
    else
        IM3FilteredSampleDelayed = IM3FilteredSample;
    end

    ErrorSignal(Sample) = IM3FilteredSampleDelayed;
 
end

% Allign the feedback error signal with the basis function signal to
% estimate the loop delay
ccorr=@(x,y) ifft(conj(fft(x)).*fft(y));
ReferenceSignal = IM3GeneratedSignal(1:length(ErrorSignal)).';
[~, LoopDelay] = max(abs(ccorr(ReferenceSignal,ErrorSignal)));
