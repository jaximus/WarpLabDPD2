function LMSFilterTaps = IM3_BlockDecorrDPD(PA_InputSignal,IM3_Basis_Orthogonal,MemoryLessPA,MemoryLessDPD, ...
    SystemFs,Signal_Bandwidth,IM3_Freq,LoopDelay,AdditionalDelay,...
    DPD_LearningBlockSize,DPD_FilteringBlockSize,nodes,RF_TX,RF_RX,node_tx,node_rx,eth_trig,Ts, Mu, NumSamples);

NumberOfBasisFunctions = length(IM3_Basis_Orthogonal(:,1));

% Plotting variables
%iq_range    = 1;                       % Plots IQ values in the range:  [-1, 1]
%rssi_range  = 1024;                    % Plots RSSI values in the range:  [0, 1024]
USE_PREAMBLE            = 1;
warp_PA_delay           = 44;          %Delat only used if USE_PREAMBLE = 0
LTS_CORR_THRESH         = 0.8;         % Normalized threshold for LTS correlation
DO_APPLY_CFO_CORRECTION = 0;           % Enable CFO estimation/correction
FFT_OFFSET              = 4;

%% Make the Third order IMD extraction filter, all frequency values are in MHz.
Fs    = round(SystemFs/1e6);    % Sampling Frequency
N     = 200;                    % Order
Fpass = (3*Signal_Bandwidth)/2; % Passband Frequency
Fstop = (5*Signal_Bandwidth)/2; % Stopband Frequency

% Calculate the coefficients using the FIRLS function.
IM3Filter  = firls(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0]);

%% Define the preamble
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

% Use 30 copies of the 16-sample STS for extra AGC settling margin
preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];

%% Adaptive Decorrelating DPD using MBF
if MemoryLessDPD
    NumBlocks = floor(NumSamples/DPD_FilteringBlockSize);
    AdaptiveFilterDelay = 1;
    LMSInput = IM3_Basis_Orthogonal;
    W = zeros(NumberOfBasisFunctions,AdaptiveFilterDelay);
    DPD_Coeff = zeros(NumBlocks,NumberOfBasisFunctions);
else
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
    payload = PA_InBlock;
    
    if(USE_PREAMBLE)
        txData  = vertcat(preamble.', payload);
    else
        txData  = payload;
    end
    
    txLength = length(txData);
    
    % Set capture lengths
    RXLength    = txLength+1000;
    
    wl_basebandCmd(nodes, 'tx_length', txLength);
    wl_basebandCmd(nodes, 'rx_length', RXLength);
    
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
    %rx_IQ = detrend(rx_IQ);
    
    % Disable the buffers and RF interfaces for TX / RX
    wl_basebandCmd(nodes, 'RF_ALL', 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, 'RF_ALL', 'tx_rx_dis');
    
    
    %% Correlate for LTS
    if(USE_PREAMBLE)
        % Complex cross correlation of Rx waveform with time-domain LTS
        lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx_IQ)));
        
        % Skip early and late samples
        lts_corr = lts_corr(32:end-32);
        
        % Find all correlation peaks
        lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));
        
        % Select best candidate correlation peak as LTS-payload boundary
        [LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
        [lts_second_peak_index,y] = find(LTS2-LTS1 == length(lts_t));
        
        % Stop if no valid correlation peak was found
        if(isempty(lts_second_peak_index))
            fprintf('No LTS Correlation Peaks Found!\n');
            return;
        end
        
        % Set the sample indices of the payload symbols and preamble
        payload_ind = lts_peaks(max(lts_second_peak_index))+32;
        lts_ind = payload_ind-160;
        
        if(DO_APPLY_CFO_CORRECTION)
            %Extract LTS (not yet CFO corrected)
            rx_lts = rx_IQ(lts_ind : lts_ind+159);
            rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
            rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);
            
            %Calculate coarse CFO est
            rx_cfo_est_lts = mean(unwrap(angle(rx_lts1 .* conj(rx_lts2))));
            rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
        else
            rx_cfo_est_lts = 0;
        end
        
        % Apply CFO correction to raw Rx waveform
        rx_cfo_corr_t = exp(1i*2*pi*rx_cfo_est_lts*[0:length(rx_IQ)-1]);
        rx_dec_cfo_corr = rx_IQ .* rx_cfo_corr_t.';
        
        % Re-extract LTS for channel estimate
        rx_lts = rx_dec_cfo_corr(lts_ind : lts_ind+159);
        rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
        rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);
        
        rx_lts1_f = fft(rx_lts1);
        rx_lts2_f = fft(rx_lts2);
        
        % Calculate channel estimate
        rx_H_est = lts_f .* (rx_lts1_f.' + rx_lts2_f.')/2;
        
        % Extract the payload samples
        payload_vec = rx_dec_cfo_corr(payload_ind : payload_ind+length(payload)-1);
    else
        payload_vec = rx_IQ(warp_PA_delay:warp_PA_delay+DPD_FilteringBlockSize-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize results FOR DEBUGGING PURPOSES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     figure(1);clf;
    %     ax(1) = subplot(2,1,1);
    %     plot(0:(length(rx_IQ)-1),real(rx_IQ))
    %     xlabel('Sample Index')
    %     title('Received I')
    %     axis([1 RXLength -iq_range iq_range])
    %
    %     ax(2) = subplot(2,1,2);
    %     plot(0:(length(rx_IQ)-1),imag(rx_IQ))
    %     xlabel('Sample Index')
    %     title('Received Q')
    %     axis([1 RXLength -iq_range iq_range])
    %
    %     linkaxes(ax,'xy')
    
    PA_OutBlock = payload_vec;
    
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
    figure();plot(TimeAxis,real(FinalCoeff(:,1)),'k');hold on;
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
    
    figure();plot(TimeAxis,real(FinalCoeff(:,1)),'k');hold on;
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


