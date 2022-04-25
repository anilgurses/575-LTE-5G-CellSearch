clear all;


addpath(genpath('./LTE'));
addpath(genpath('./5G'));

% Initializing simulation parameters and modes
Modulations = ["QPSK" "16QAM" "64QAM" "256QAM"];
ChModels = ["awgn" "rayleigh"];
modes = ["5G"];
nrSSbursts = ["Case A" "Case B" "Case C" "Case D" "Case E"];
nrminBws = [5 10 40];
lteRbs = [6 15 25 50 100];


for snr = -20:10:40
    for mod = 1:length(modes)
        for ch = 1:length(ChModels)
            for m = 1:length(Modulations)
                % Getting configurations
                modul = Modulations(m);
                mode = modes(mod);
                chan = ChModels(ch);

                if mode == "LTE"
                    rmc = lteRMCDL('R.3');
                    rmc.PDSCH.Modulation{1} = modul;
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
                    lteCellSearchSync(6,rmc, 'awgn', [], [], 30, false);
                else
                    for nss = 1:numel(nrSSbursts)
                        for bw = 1:numel(nrminBws)
                            nrSS = nrSSbursts(nss);
                            bandw = nrminBws(bw);
                            fnRoot = sprintf('Results/NR/%s/%s/%s/BW-%dMHz/SNR(%ddB)', chan, modul, nrSS, bandw, snr);
                            if ~exist(fnRoot, 'dir')
                              mkdir(fnRoot);
                            end
                            fn = sprintf("%s/logs.txt", fnRoot);
                            diary(fn);
                            diary("on");
                            
                            lg = sprintf("Simulation started. Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                            disp(lg);
                            config = struct();
                            config.NCellID = 102;
                            
                            % Configure an SS burst
                            config.BlockPattern = nrSS;         % FR1: 'Case A','Case B','Case C'. FR2: 'Case D','Case E'
                            config.TransmittedBlocks = ones(1,8);   % Bitmap of SS blocks transmitted
                            config.SubcarrierSpacingCommon = 15;    % SIB1 subcarrier spacing in kHz (15 or 30 for FR1. 60 or 120 for FR2)
                            config.EnableSIB1 = 1;                  % Set to 0 to disable SIB1
                            config.Modulation = modul;
                            
                            % Set the minimum channel bandwidth for the NR band required to
                            % configure CORESET0 in FR1 (See TS 38.101-1 Table 5.3.5-1)
                            config.MinChannelBW = bandw; % 5, 10, 40 MHz
                        
                            ssbIdx = 0; % Index of the SSB to boost (0-based)
                            boost = 6; % SNR boost in dB
                            nrbSSB = 20;
                            SNRdB = snr; % SNR for AWGN
                        
                            
                            try
                                wavegenConfig = hSIB1WaveformConfiguration(config);
                                NRCellSearchSync(wavegenConfig, config, 'awgn', [], [], ssbIdx, boost, nrbSSB, SNRdB, false);
                            catch
                                figHandles = findall(0,'Type','figure');
                                err = sprintf("ERROR --- Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                                warning(err);
                                disp("Encountered an error");
                                close(figHandles);
                                diary("off"); 
                            end

                            figHandles = findall(0,'Type','figure'); 
                            
                             for i = 1:numel(figHandles)
                                 pth = sprintf('%s/fig-%d.png', fnRoot, i);
                                 saveas(figHandles(i),pth,'png')
                                 pthm = sprintf('%s/fig-%d.m', fnRoot, i);
                                 saveas(figHandles(i),pthm,'m')
                             end
                             close(figHandles);
                             diary("off"); 
                        end
                    end
                end
            end
        end
    end
end
