clear all;
close all;
clc;

% rand('state',0);
% set(0,'defaultlinelinewidth',1.5)

%% Simulation Paramters %%
TxScenario                          =   3;          % 1- Contigious Single CC
% 2- Multi-Cluster Single CC
% 3- Intra-Band CA Two CC
N_symbols                           =  2000;        % No of OFDM Symbols for Simulation
TxSignalType                        =   2;          % 1- OFDM for DL Transmission
% 2- SCFDMA for UL Transmission
ModulationType                      =   1;          % 1- QPSK
% 2- 16-QAM
% 3- 64-QAM
POWER_PLOT_1MHZ                     =   0;          % Use 1 MHz power plot

% RF and PA paramters
Pout_dBm                            =   21;         % Desired Tx Output Power after the PA in dBm

MemoryLessPA                        =   1;          % 1- MemoryLess PA model
% 0- PH PA model with memory
MemoryLessDPD                       =   1;          % 1- MemoryLess DPD
% 0- MemoryDPD

USE_WARP = 0;                                       %For a simulation only. 
DoTraining = 1;                                     %To set alpha yourself. No training
IM3_BlockDecorrDPD_Coeffs           = 0.5+0.1*i;    %For Memoryless Pretraining

DPD_LearningBlockSize  = 500;                       % Decorrelating DPD Learning Block size
DPD_FilteringBlockSize = 1000;                      % Decorrelating DPD Filtering Block size
NumSamples = 2000000;                               % Total number of samples used for learning
Mu = 4;                                             %LMS Gain

ScalingForPA = 4.2;                                 %Changes magnitude of PA input signal

%% Power Amplifier Model
PA_Power_Measured = 23; % The Tx power at which the PA paramters were measured at
% PH PA Branch Filters of the used PH PA model
PH_f1 = [0.9741 + 0.1014i;
    0.2833 - 0.3192i;
    -0.5420 + 0.1199i;
    0.6692 + 0.3308i;
    -0.4427 - 0.4366i;
    0.1002 + 0.1407i];

PH_f3 = [-0.0489 + 0.0192i;
    0.0025 + 0.0063i;
    -0.0093 - 0.0098i;
    0.0080 + 0.0049i;
    -0.0048 - 0.0028i;
    0.0024 - 0.0006i];

PH_f5 = [0.0015 - 0.0029i;
    0.0000 + 0.0002i;
    0.0000 + 0.0002i;
    0.0000 + 0.0003i;
    0.0000 + 0.0002i;
    0.0000 + 0.0002i];

% The memoryless paramters of the PA model (Not used in this m-file)
MemorylessPA_Paramters = [1.0735 - 0.0287i;
    0.5 + 0.5;
    0.0012 - 0.0030i];
Beta_1 = MemorylessPA_Paramters(1);
Beta_3 = MemorylessPA_Paramters(2);
Beta_5 = MemorylessPA_Paramters(3);


%% Baseband Tx signal paramters
LTE_Bandwidth = 1.4; % Carrier BW of the LTE Tx Signal
NRB1 = 6;           % Number of RB allocated to first CC
NRB2 = 6;           % Number of RB allocated to second CC
CarrierSpacing = 6; % Spacing in MHz between the 2 CC
IM3_Freq = 3*(CarrierSpacing/2);
Signal_Bandwidth = NRB1*0.18;

%% Setup WARPLAB
if(USE_WARP)
    NUMNODES       = 1;
    % TX variables
    MAX_TX_LEN     = 2^25;                 % 2^14 =    16384 --> Max TX / RX length for WARP v2
    % 2^15 =    32768 --> Max TX / RX length for WARP v3 (WARPLab 7.4.0 and prior)
    % 2^20 =  1048576 --> Soft max TX / RX length for WARP v3 Java Transport (WARPLab 7.5.x)
    % 2^25 = 33554432 --> Soft max TX / RX length for WARP v3 MEX Transport (WARPLab 7.5.x)
    % RX variables
    USE_AGC        = false;
    ManualRxGainRF = 2;                    % Rx RF Gain in [1:3] (ignored if USE_AGC is true)
    ManualRxGainBB = 9;                   % Rx Baseband Gain in [0:31] (ignored if USE_AGC is true)
    
    % TX variables
    BB_GAIN = 3;                           %Must be integer in [0,1,2,3] for approx ![-5, -3, -1.5, 0]dB baseband gain
    RF_GAIN = 41;                          %Must be integer in [0:63] for approx [0:31]dB RF gain
    
    % Create a vector of node objects
    nodes = wl_initNodes(NUMNODES);
    
    % If we only have one node, then we should do RF loopback.
    % Otherwise, we should transmit from nodes(1) to nodes(2)
    if (NUMNODES == 1)
        node_tx = nodes(1);
        node_rx = nodes(1);
    else
        node_tx = nodes(1);
        node_rx = nodes(2);
    end
    
    %S Create a UDP broadcast trigger and tell each node to be ready for it
    eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes, 'add_ethernet_trigger', [eth_trig]);
    
    % Read Trigger IDs into workspace
    [T_IN_ETH_A, T_IN_ENERGY, T_IN_AGCDONE, T_IN_REG, T_IN_D0, T_IN_D1, T_IN_D2, T_IN_D3, T_IN_ETH_B] =  wl_getTriggerInputIDs(nodes(1));
    [T_OUT_BASEBAND, T_OUT_AGC, T_OUT_D0, T_OUT_D1, T_OUT_D2, T_OUT_D3] = wl_getTriggerOutputIDs(nodes(1));
    
    % For both nodes, we will allow Ethernet to trigger the buffer baseband and the AGC
    wl_triggerManagerCmd(nodes, 'output_config_input_selection', [T_OUT_BASEBAND, T_OUT_AGC], [T_IN_ETH_A]);
    
    % Set the trigger output delays.
    %
    % NOTE:  We are waiting 3000 ns before starting the AGC so that there is time for the inputs
    %   to settle before sampling the waveform to calculate the RX gains.  If you have different
    %   timing between nodes, then you will need to adjust this timing so that you sample the
    %   waveform in the correct spot to set valid RX gains.  For example, if Node 1 was a WARP v2
    %   node and Node 1 was a WARP v3 node, then we will not be able to use the standard 3000 ns
    %   since it takes WARP v2 using software Ethernet triggers about 25 us longer than WARP v3
    %   using hardware Ethernet triggers (this was observed experimentally using fixed gain values).
    %   Therefore, in this case, we need to add 25000 ns to the standard 3000 ns to give an AGC
    %   delay of 28000 ns so that we sample the waveform in the correct spot to compute valid
    %   RX gains.
    %
    nodes.wl_triggerManagerCmd('output_config_delay', [T_OUT_BASEBAND], 0);
    nodes.wl_triggerManagerCmd('output_config_delay', [T_OUT_AGC], 3000);     % 3000 ns delay before starting the AGC
    
    
    % Get IDs for the interfaces on the boards. Since this example assumes each
    % board has the same interface capabilities, we only need to get the IDs
    % from one of the boards
    [RFA,RFB] = wl_getInterfaceIDs(nodes(1));
    
    % If we only have one node, then we should do RF loopback from RFA to RFB.
    % Otherwise, we should use the RFA interface of each node.
    if (NUMNODES == 1)
        RF_TX = RFB;
        RF_RX = RFA;
    else
        RF_TX = RFA;
        RF_RX = RFA;
    end
    
    % Set up the interface for the experiment
    wl_interfaceCmd(nodes, 'RF_ALL', 'tx_gains', BB_GAIN, RF_GAIN);
    wl_interfaceCmd(nodes, 'RF_ALL', 'channel', 2.4, 6);
    
    if(USE_AGC)
        wl_interfaceCmd(nodes, 'RF_ALL', 'rx_gain_mode', 'automatic');
        wl_basebandCmd(nodes, 'agc_target', -10);
    else
        wl_interfaceCmd(nodes, 'RF_ALL', 'rx_gain_mode', 'manual');
        wl_interfaceCmd(nodes, 'RF_ALL', 'rx_gains', ManualRxGainRF, ManualRxGainBB);
    end
    
    %First, read the transmitter's maximum I/Q buffer length.  This example assumes
    % that each board has the same interface capabilities, so we only need to read
    % the max I/Q buffer length of node_tx RFA.
    maximum_buffer_len = wl_basebandCmd(node_tx, RF_TX, 'tx_buff_max_num_samples');
    
    Ts = 1/(wl_basebandCmd(nodes(1),'tx_buff_clk_freq'));
    Ts_RSSI = 1/(wl_basebandCmd(nodes(1),'rx_rssi_clk_freq'));
end

%% Baseband Equivilent LTE Signal Transmitter
[LTE_Signal CC1 CC2 SystemFs UpsamplingFactor] = LTE_Transmitter(LTE_Bandwidth,CarrierSpacing,NRB1,NRB2,...
    N_symbols,TxSignalType,ModulationType,TxScenario);
% Scale the Baseband generated signal to have a unit RMS power
% TX_PowerScale = sqrt(10^((PA_Power_Measured-Pout_dBm)/10));
% RMS_PAin = sqrt(mean(abs(LTE_Signal).^2));
% ScalingForPA = 1/(RMS_PAin*TX_PowerScale);

PA_InputSignal = LTE_Signal*ScalingForPA;
CC1 = CC1*ScalingForPA;
CC2 = CC2*ScalingForPA;

%% Broadcast without DPD
if(USE_WARP)
    txData  = [PA_InputSignal];
    txLength = length(txData);
    
    % Set capture lengths
    RXLength    = txLength;
    wl_basebandCmd(nodes, 'tx_delay', 0);
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
    
    % Disable the buffers and RF interfaces for TX / RX
    wl_basebandCmd(nodes, 'RF_ALL', 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, 'RF_ALL', 'tx_rx_dis');
    PA_OutputSignal = rx_IQ(50:end);         %throw away 50 early samples
else
    nodes = 0;RF_TX = 0;RF_RX = 0;,node_tx = 0;,node_rx = 0;,eth_trig = 0;,Ts = 0;
    % Power Amplifier Model
    if MemoryLessPA
        PA_OutputSignal = MemoryLess_PA(PA_InputSignal,Beta_1,Beta_3,Beta_5);
    else
        PA_OutputSignal = MemoryPH_PA(PA_InputSignal,PH_f1,PH_f3,PH_f5);
    end
end
%% Decorrelating DPD
% Generate the IM3 third order basis function
IM3GeneratedSignal = (conj(CC2).*(CC1.^2));
% Additional delay added to test the DPD performance with any delay
% that can exist in real implementation
AdditionalDelay_Musecs = 0; % In MicroSeconds
AdditionalDelay_Samples = round(AdditionalDelay_Musecs*1e-6*SystemFs);

%% Loop Delay Estimation
LoopDelay = DPD_LoopDelayEst(PA_InputSignal, IM3GeneratedSignal, MemoryLessPA, SystemFs,Signal_Bandwidth, ...
    IM3_Freq,Beta_1,Beta_3,Beta_5,PH_f1,PH_f3,PH_f5,AdditionalDelay_Samples);

%% IM3 Block Decorrelating DPD Estimation
if (DoTraining == 1)
    IM3_BlockDecorrDPD_Coeffs = IM3_BlockDecorrDPD(PA_InputSignal,IM3GeneratedSignal.',MemoryLessPA,MemoryLessDPD, ...
        SystemFs,Signal_Bandwidth,IM3_Freq, ...
        LoopDelay,AdditionalDelay_Samples,DPD_LearningBlockSize,DPD_FilteringBlockSize,nodes,RF_TX,RF_RX,node_tx,node_rx,eth_trig,Ts,Mu,NumSamples,USE_WARP,...
        Beta_1,Beta_3,Beta_5,PH_f1,PH_f3,PH_f5)
end
%% Applying IM3 Decorrelating DPD
AdaptiveFilterDelay  = length(IM3_BlockDecorrDPD_Coeffs) - 1;
IM3_DPD_Signal = conv(IM3GeneratedSignal,IM3_BlockDecorrDPD_Coeffs);
IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
n = (1:length(IM3_DPD_Signal)).';
IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
PAin_IM3_DPD   = PA_InputSignal + IM3_DPD_Signal;

%% Broadcast with DPD
if(USE_WARP)
    txData  = [PAin_IM3_DPD];
    txLength = length(txData);
    
    % Set capture lengths
    RXLength    = txLength;
    
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
    
    % Disable the buffers and RF interfaces for TX / RX
    wl_basebandCmd(nodes, 'RF_ALL', 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, 'RF_ALL', 'tx_rx_dis');
    
    PA_Output_IM3_DPD = rx_IQ(50:end);         %throw away 50 early samples
else
    % Power Amplifier Model
    if MemoryLessPA
        PA_Output_IM3_DPD = MemoryLess_PA(PAin_IM3_DPD,Beta_1,Beta_3,Beta_5);
    else
        PA_Output_IM3_DPD = MemoryPH_PA(PAin_IM3_DPD,PH_f1,PH_f3,PH_f5);
    end
end
%% Spectral plots
figure(3);
plot_freqdomain(PA_OutputSignal,SystemFs,'','r',UpsamplingFactor,POWER_PLOT_1MHZ,23);
hold on;
plot_freqdomain(PA_Output_IM3_DPD,SystemFs,'','k',UpsamplingFactor,POWER_PLOT_1MHZ,23);
grid on;
L2= legend('No DPD','$IM3$ Decorr.');
set(L2,'Interpreter','Latex');
