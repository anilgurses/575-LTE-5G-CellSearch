addpath('./LTE');
addpath('./5G');


mode = "LTE";

if mode == "LTE"
    rmc = lteRMCDL('R.3');
    rmc.NCellID = 17;
    rmc.TotSubframes = 41; % Resource 
    rmc.PDSCH.RNTI = 61;
    % SIB parameters
    rmc.SIB.Enable = 'On';
    rmc.SIB.DCIFormat = 'Format1A';
    rmc.SIB.AllocationType = 0;
    rmc.SIB.VRBStart = 8;
    rmc.SIB.VRBLength = 8;
    rmc.SIB.Gap = 0;
    % SIB data field filled with random bits, this is not a valid SIB
    % message
    rmc.SIB.Data = randi([0 1],176,1);
    lteCellSearchSync(6,rmc, 'awgn', 30);

else
    config = struct();
    config.NCellID = 102;
    
    % Configure an SS burst
    config.BlockPattern = 'Case B';         % FR1: 'Case A','Case B','Case C'. FR2: 'Case D','Case E'
    config.TransmittedBlocks = ones(1,8);   % Bitmap of SS blocks transmitted
    config.SubcarrierSpacingCommon = 15;    % SIB1 subcarrier spacing in kHz (15 or 30 for FR1. 60 or 120 for FR2)
    config.EnableSIB1 = 1;                  % Set to 0 to disable SIB1
    
    % Set the minimum channel bandwidth for the NR band required to
    % configure CORESET0 in FR1 (See TS 38.101-1 Table 5.3.5-1)
    config.MinChannelBW = 5; % 5, 10, 40 MHz

    ssbIdx = 0; % Index of the SSB to boost (0-based)
    boost = 6; % SNR boost in dB
    nrbSSB = 20;
    SNRdB = 20; % SNR for AWGN

    NRCellSearchSync(config,ssbIdx, boost, nrbSSB, SNRdB);
end