clear all;


addpath(genpath('./LTE'));
addpath(genpath('./5G'));

% Initializing simulation parameters and modes
Modulations = ["QPSK" "16QAM" "64QAM" "256QAM"];
ChModels = ["awgn" "rayleigh"];
modes = ["5G" "LTE"];
nrSSbursts = ["Case A" "Case B" "Case C" "Case D" "Case E"];
nrMaxBws = [25 100 150 200 300 400];
nrRbs = [50 100 125 200 275];
lteRbs = [6 15 25 50 100];
numberOfPorts = [1 2 4];


for snr = -10:10:40
    for mod = 1:length(modes)
        for ch = 1:length(ChModels)
            for m = 1:length(Modulations)
                % Getting configurations
                modul = Modulations(m);
                mode = modes(mod);
                chan = ChModels(ch);
                if (mode == "LTE")
                    for lrb = 1:numel(lteRbs)
                        dlrb = lteRbs(lrb);
                        bw = "1_4MHz";
                        switch dlrb
                            case 6
                                bw = "1_4MHz";
                            case 15
                                bw = "3MHz";
                            case 25
                                bw = "5MHz";
                            case 50
                                bw = "10MHz";
                            case 100
                                bw = "20MHz";
                        end

                        fnRoot = sprintf('Results/LTE/%s/%s/%s/SNR(%ddB)', chan, modul, bw, snr);

                        if ~exist(fnRoot, 'dir')
                          mkdir(fnRoot);
                        end
                        fn = sprintf("%s/logs.txt", fnRoot);
                        diary(fn);
                        diary("on");
                        
                        lg = sprintf("LTE Cell Search Simulation started. Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                        disp(lg);
                        
                        rmc = lteRMCDL("R.0");
                        rmc.PDSCH.Modulation{1} = modul;
                        rmc.NCellID = 17;
                        rmc.TotSubframes = 41; 
                        rmc.PDSCH.RNTI = 61;
    
                        rmc.NDLRB = dlrb; % Resource block  
                        % SIB parameters
                        rmc.SIB.Enable = 'On';
                        rmc.SIB.DCIFormat = 'Format1A';
                        rmc.SIB.AllocationType = 0;
                        rmc.SIB.VRBStart = 2;
                        rmc.SIB.VRBLength = 2;
                        rmc.SIB.Gap = 0;
                        % SIB data field filled with random bits, this is not a valid SIB
                        % message
                        rmc.SIB.Data = randi([0 1],176,1);

                        rmcoverrided = lteRMCDL(rmc);

                        try
                            lteCellSearchSync(rmcoverrided, chan, [0 0.00015], [2 3], snr, fnRoot, false);
                        catch e
                            figHandles = findall(0,'Type','figure');
                            err = sprintf("ERROR --- Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                            warning(err);
                            warning(e.message);
                            disp("Encountered an error");
                            disp(e.message);
                            close(figHandles);
                            diary("off"); 
                        end

                        figHandles = findall(0,'Type','figure'); 
                            
                         for i = 1:numel(figHandles)
                             pth = sprintf('%s/fig-%d.png', fnRoot, i);
                             saveas(figHandles(i),pth,'png')
                             pthm = sprintf('%s/fig-%d.fig', fnRoot, i);
                             saveas(figHandles(i),pthm,'fig')
                         end

                         close(figHandles);
                         diary("off"); 
                    end
                else
                    for nss = 1:numel(nrSSbursts)
                        for rb = 1:numel(nrRbs)
                            nrSS = nrSSbursts(nss);
                            resBlk = nrRbs(rb);

                            config = struct();
                            config.NCellID = 102;
                            
                            % Configure an SS burst
                            config.BlockPattern = nrSS;         % FR1: 'Case A','Case B','Case C'. FR2: 'Case D','Case E'
                            config.TransmittedBlocks = ones(1,8);   % Bitmap of SS blocks transmitted
                            config.SubcarrierSpacingCommon = 15;    % SIB1 subcarrier spacing in kHz (15 or 30 for FR1. 60 or 120 for FR2)
                            config.EnableSIB1 = 1;                  % Set to 0 to disable SIB1
                            config.Modulation = modul;
                            ssbSpace = hSSBurstSubcarrierSpacing(config.BlockPattern);

                            % Number of resource blocks
                          
                            if (ssbSpace <= 240 && resBlk > 131)
                                continue
                            end
                            config.NSizeGrid = resBlk; 

                            bw = 0;
                                                        
                            % Set the minimum channel bandwidth for the NR band required to
                            % configure CORESET0 in FR1 (See TS 38.101-1 Table 5.3.5-1)
                            switch ssbSpace
                                case 15
                                    config.MinChannelBW = 5;
                                    config.Bandwidth = 50;
                                    bw = 0.18*resBlk;
                                case 30
                                    config.MinChannelBW = 10;
                                    config.Bandwidth = 100;
                                    bw = 0.36*resBlk;
                                case 60
                                    config.MinChannelBW = 20;
                                    config.Bandwidth = 200;
                                    bw = 0.72*resBlk;
                                case 120
                                    config.MinChannelBW = 40;
                                    config.Bandwidth = 400;
                                    bw = 1.44*resBlk;
                                case 240
                                    config.MinChannelBW = 70;
                                    config.Bandwidth = 400;
                                    bw = 2.88*resBlk;
                            end
                           
                            fnRoot = sprintf('Results/NR/%s/%s/%s/%dMHz/SNR(%ddB)', chan, modul, nrSS, bw, snr);

                            if ~exist(fnRoot, 'dir')
                              mkdir(fnRoot);
                            else
                                continue
                            end
                            fn = sprintf("%s/logs.txt", fnRoot);
                            diary(fn);
                            diary("on");
                            
                            lg = sprintf("5G NR Cell Search Simulation started. Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                            disp(lg);
                            
                            ssbIdx = 0; % Index of the SSB to boost (0-based)
                            boost = 6; % SNR boost in dB
                            nrbSSB = 20;
                            SNRdB = snr; % SNR for AWGN
                           
                            try
                                wavegenConfig = hSIB1WaveformConfiguration(config);
                                NRCellSearchSync(wavegenConfig, config, chan, [0 0.00015], [2 3], ssbIdx, boost, nrbSSB, SNRdB, fnRoot, false);
                            catch e
                                figHandles = findall(0,'Type','figure');
                                err = sprintf("ERROR --- Params -> Channel: %s | Modulation: %s | SNR: %d",chan, modul, snr);
                                warning(err);
                                warning(e.message);
                                disp(e.message);
                                disp("Encountered an error");
                                close(figHandles);
                                diary("off"); 
                            end

                            figHandles = findall(0,'Type','figure'); 
                            
                             for i = 1:numel(figHandles)
                                 pth = sprintf('%s/fig-%d.png', fnRoot, i);
                                 saveas(figHandles(i),pth,'png')
                                 pthm = sprintf('%s/fig-%d.fig', fnRoot, i);
                                 saveas(figHandles(i),pthm,'fig')
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
