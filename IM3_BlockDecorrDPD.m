function LMSFilterTaps = IM3_BlockDecorrDPD(PA_InputSignal,IM3_Basis_Orthogonal,MemoryLessPA,MemoryLessDPD, ...
    SystemFs,Signal_Bandwidth,IM3_Freq,LoopDelay,AdditionalDelay,...
    DPD_LearningBlockSize,DPD_FilteringBlockSize,nodes,RF_TX,RF_RX,node_tx,node_rx,eth_trig,Ts);

NumberOfBasisFunctions = length(IM3_Basis_Orthogonal(:,1));

% Third order IMD extraction filter, all frequency values are in MHz.
Fs    = round(SystemFs/1e6);  % Sampling Frequency
N     = 100;  % Order
Fpass = (3*Signal_Bandwidth)/2; % Passband Frequency
Fstop = (5*Signal_Bandwidth)/2; % Stopband Frequency
Wpass = 1;    % Passband Weight
Wstop = 5;    % Stopband Weight
dens  = 20;   % Density Factor
% Calculate the coefficients using the FIRPM function.
IM3Filter  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
    {dens});

% Adaptive Decorrelating DPD using MBF
if MemoryLessDPD
    NumSamples = 1000000; % Total number of samples used for learning
    NumBlocks = floor(NumSamples/DPD_FilteringBlockSize);
    Mu = 1;
    AdaptiveFilterDelay = 1;
    LMSInput = IM3_Basis_Orthogonal;
    W = zeros(NumberOfBasisFunctions,AdaptiveFilterDelay);
    DPD_Coeff = zeros(NumBlocks,NumberOfBasisFunctions);
else
    NumSamples = 20000000;
    Mu = 0.5;
    NumBlocks = floor(NumSamples/DPD_FilteringBlockSize);
    AdaptiveFilterDelay = 5;
    LMSInput = [IM3_Basis_Orthogonal, zeros(1,AdaptiveFilterDelay-1)];
    W = zeros(1,AdaptiveFilterDelay);
    DPD_Coeff = zeros(NumBlocks,AdaptiveFilterDelay);
    Alpha = zeros(AdaptiveFilterDelay,1);
end

% Reset the delay line
DelayLineReset = 1;
DelayLine(0,0,DelayLineReset);
DelayLineReset = 0;

DPD_BlockIndx = 0;

for Sample = 1:DPD_FilteringBlockSize:NumSamples
    
    PA_In_StartIndx = Sample;
    PA_In_EndIndx   = Sample + DPD_FilteringBlockSize - 1;
    DPD_BlockIndx   = DPD_BlockIndx + 1;
    
    if Sample > (LoopDelay + DPD_FilteringBlockSize) + AdaptiveFilterDelay - 1
        
        if MemoryLessDPD
            
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = Sample - (LoopDelay + DPD_FilteringBlockSize) + 1;
            LMS_In_EndIndx = LMS_In_StartIndx + DPD_LearningBlockSize - 1;
            LMS_InputBlock = LMSInput(:,LMS_In_StartIndx:LMS_In_EndIndx);
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/DPD_LearningBlockSize;
            
            % Block LMS filter update
            W = W - Mu.*MeanCorrelation;
            
            % Store DPD filter coeff.
            DPD_Coeff(DPD_BlockIndx,:) = W';
            
            % Apply LMS block filtering
            Alpha1 = W';
            LMS_FilteringBlock = LMSInput(:,PA_In_StartIndx:PA_In_EndIndx);
            IM3_Block  = (Alpha1.*LMS_FilteringBlock).';
            
            % Add the LMS filter output to the PA input
            PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx) ...
                + IM3_Block.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*IM3_Freq*1e6/SystemFs);
            
        else
            
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = Sample - (LoopDelay + DPD_FilteringBlockSize) + 1;
            LMS_In_EndIndx = LMS_In_StartIndx + DPD_LearningBlockSize - 1;
            
            LMS_InputBlock = [LMSInput(LMS_In_StartIndx:LMS_In_EndIndx); ...
                LMSInput(LMS_In_StartIndx - 1:LMS_In_EndIndx - 1); ...
                LMSInput(LMS_In_StartIndx - 2:LMS_In_EndIndx - 2); ...
                LMSInput(LMS_In_StartIndx - 3:LMS_In_EndIndx - 3); ...
                LMSInput(LMS_In_StartIndx - 4:LMS_In_EndIndx - 4)];
            
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/DPD_LearningBlockSize;
            
            IM3_Block = zeros(DPD_FilteringBlockSize,1);
            for Tap = 1:AdaptiveFilterDelay
                % Block LMS filter update
                W(Tap) = W(Tap) - Mu.*MeanCorrelation(Tap);
                % Store DPD filter coeff.
                DPD_Coeff(DPD_BlockIndx,Tap) = W(Tap)';
                % Apply LMS block filtering
                Alpha(Tap) = W(Tap)';
                LMS_FilteringBlock = LMSInput(PA_In_StartIndx-(Tap-1):PA_In_EndIndx-(Tap-1));
                IM3_Block = IM3_Block + Alpha(Tap).*LMS_FilteringBlock.';
            end
            
            % Add the LMS filter output to the PA input
            PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx) ...
                + IM3_Block.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*IM3_Freq*1e6/SystemFs);
            
        end
    else
        PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx);
    end
    
    % Broadcast Block using WARPLAB
    PA_InBlock;
    txData  = [PA_InBlock];
    txLength = length(txData);
    
    % Set capture lengths
    RXLength    = txLength+4000;
    
    wl_basebandCmd(nodes, 'tx_delay', 0);
    wl_basebandCmd(nodes, 'tx_length', txLength);
    wl_basebandCmd(nodes, 'rx_length', txLength);
    
    % Transmit IQ data to the TX node
    wl_basebandCmd(node_tx, [RF_TX], 'write_IQ', txData(:));
    
    % Enabled the RF interfaces for TX / RX
    wl_interfaceCmd(node_tx, RF_TX, 'tx_en');
    wl_interfaceCmd(node_rx, RF_RX, 'rx_en');
    
    % Enable the buffers for TX / RX
    wl_basebandCmd(node_tx, RF_TX, 'tx_buff_en');
    wl_basebandCmd(node_rx, RF_RX, 'rx_buff_en');
    
    % Send the Ethernet trigger to start the TX / RX
    eth_trig.send();
    
    % Wait until the TX / RX is done
    pause(1.2 * txLength * Ts);
    
    % Read the IQ and RSSI data from the RX node
    rx_IQ    = wl_basebandCmd(node_rx, [RF_RX], 'read_IQ', 0, RXLength);
    rx_gains = wl_basebandCmd(node_rx, [RF_RX], 'agc_state');
    
    % Disable the buffers and RF interfaces for TX / RX
    wl_basebandCmd(nodes, 'RF_ALL', 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, 'RF_ALL', 'tx_rx_dis');
    
    PA_OutBlock = rx_IQ(43:43+DPD_FilteringBlockSize-1);
    
    % Shift the PA output such that the IM3 frequency is at baseband
    PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*IM3_Freq*1e6/SystemFs);
    
    % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
    if Sample == 1
        [IM3FilteredBlock Current_State] = filter(IM3Filter,1,PA_OutBlockShifted);
    else
        [IM3FilteredBlock Current_State] = ...
            filter(IM3Filter,1,PA_OutBlockShifted,Previous_State);
    end
    Previous_State = Current_State;
    
    % Implement any possible additional delay due to hardware latency
    if AdditionalDelay ~= 0
        IM3FilteredBlockDelayed = zeros(size(IM3FilteredBlock));
        for BlockSample = 1:length(IM3FilteredBlock)
            IM3FilteredSample = IM3FilteredBlock(BlockSample);
            IM3FilteredSampleDelayed = DelayLine(IM3FilteredSample,AdditionalDelay,DelayLineReset);
            IM3FilteredBlockDelayed(BlockSample) = IM3FilteredSampleDelayed;
        end
    else
        IM3FilteredBlockDelayed = IM3FilteredBlock;
    end
    
    % Extract the IM3 block from the PA output in the feedback receiver
    ErrorBlock = IM3FilteredBlockDelayed(1:DPD_LearningBlockSize);
end

% Plot DPD filter taps convergence
if MemoryLessDPD
    FinalCoeff = DPD_Coeff(1:NumBlocks,:);
    Samples = 1:length(FinalCoeff(:,1));
    TimeAxis = (Samples*DPD_FilteringBlockSize/SystemFs)*1e3;
    figure(1);plot(TimeAxis,real(FinalCoeff(:,1)),'k');hold on;
    plot(TimeAxis,imag(FinalCoeff(:,1)),'r')
    xlabel('Time in msecs'); ylabel('DPD Filter Tap')
    grid on;
    L2= legend('Real($\alpha_1$)','Imag($\alpha_1$)');
    set(L2,'Interpreter','Latex');
    LMSFilterTaps = Alpha1;
    
else
    
    FinalCoeff = DPD_Coeff(1:NumBlocks,:);
    Blocks = 1:length(FinalCoeff(:,1));
    TimeAxis = (Blocks*DPD_FilteringBlockSize/SystemFs)*1e3;
    
    figure(1);plot(TimeAxis,real(FinalCoeff(:,1)),'k');hold on;
    plot(TimeAxis,real(FinalCoeff(:,2)),'r');hold on;
    plot(TimeAxis,real(FinalCoeff(:,3)),'b');hold on;
    plot(TimeAxis,real(FinalCoeff(:,4)),'g');hold on;
    plot(TimeAxis,real(FinalCoeff(:,5)),'m');
    xlabel('Time in msecs'); ylabel('DPD Filter Taps')
    grid on;
    L2= legend('Real($\alpha_1(0)$)','Real($\alpha_1(1)$)','Real($\alpha_1(2)$)',...
        'Real($\alpha_1(3)$)','Real($\alpha_1(4)$)');
    set(L2,'Interpreter','Latex');
    
    LMSFilterTaps = Alpha.';
    
end


