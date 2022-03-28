%hSIB1WaveformConfiguration MIB/SIB1 configuration of nrWaveformGenerator
%   WAVEGENCONFIG = hSIB1WaveformConfiguration(CONFIG) creates an
%   nrDLCarrierConfig configuration object WAVEGENCONFIG that can be used
%   to generate an SS burst carrying a master information block and the
%   control and data channels carrying the first system information block
%   (SIB1) using the configuration CONFIG.

%   Copyright 2020-2021 The MathWorks, Inc.

function wavegenConfig = hSIB1WaveformConfiguration(config)

    % Map configuration to nrDLCarrierConfig
    wavegenConfig = nrDLCarrierConfig;
    wavegenConfig.NCellID = config.NCellID;
    
    % Configure SSBurst and MIB content related to SIB1
    ssBurst = getWavegenSSBurstConfig(config);
    transmittedBlocks = ssBurst.TransmittedBlocks;
    
    % Configure common SCS carrier and BWP for control and data using
    % the MIB content in SS burst. The BWP spans exactly CORESET 0. The
    % frequency offset between CORESET 0 and SS burst is determined by
    % the parameter PDCCH-ConfigSIB1 in MIB, as defined in TS 38.213
    % Tables 13-1 through 13-10. CORESET 0 is located in a multiple of 6
    % RB from CRB0.
    [commonCarrier,commonBwp] = getCommonSCSCarrierAndBWPConfig(ssBurst,config.MinChannelBW);
    
    % For cases of multiple numerologies, create an additional SCS
    % carrier to associate to the SS burst
    carriers{1} = commonCarrier;
    scsSSB = hSSBurstSubcarrierSpacing(ssBurst.BlockPattern);
    scsCommon = ssBurst.SubcarrierSpacingCommon;
    if scsSSB ~= scsCommon
        % Subcarrier spacing carriers and bandwidth part configuration for SSB
        ssbCarrier = nrSCSCarrierConfig;
        ssbCarrier.NSizeGrid = commonCarrier.NSizeGrid*scsCommon/scsSSB;
        ssbCarrier.SubcarrierSpacing = scsSSB;
        carriers{2} = ssbCarrier;
    end
    
    bwps{1} = commonBwp;
    
    % Configure CORESET 0 and Search Spaces
    msbIdx = floor(ssBurst.PDCCHConfigSIB1/16);
    scsPair = [scsSSB scsCommon];
    kssb = ssBurst.KSSB;
    % Read CORESET 0 configuration from TS 38.213 Tables 13-1 to 13-10
    [csetNRB,csetDuration,~,csetPattern] = hCORESET0Resources(msbIdx,scsPair,config.MinChannelBW,kssb);
    
    coreset = nrCORESETConfig();
    coreset.FrequencyResources = ones(1,csetNRB/6);
    coreset.Duration = csetDuration;
    
    % Set other parameters of CORESET 0 according to TS 38.211 Section 7.3.2.2
    coreset.CCEREGMapping = 'interleaved';
    coreset.REGBundleSize = 6;
    coreset.InterleaverSize = 2;
    coreset.ShiftIndex = config.NCellID;
    coreset.CORESETID = 0;
    
    % Configure as many search spaces as SS blocks in burst
    ss = getType0CommonSearchSpaces(ssBurst,csetNRB,csetDuration,csetPattern);
    
    % Construct DCI Format 1_0 message for PDSCH generation
    % Default configuration
    dci = DCIFormat1_0_SIRNTI(csetNRB);
    dci.FrequencyDomainResources = hRIV(csetNRB,0,8);   % RIV: NSizeBWP, RBStart and LRBS associated to PDSCH
    dci.TimeDomainResources = 0;                        % Time domain resource assignment (row index - 1) in TS 38.214 Table 5.1.2.1.1-2
    dci.VRBToPRBMapping = 0;                            % Virtual resource blocks to physical resource blocks interleaving (1/0)
    dci.ModulationCoding = 0;                           % Modulation and coding scheme (0-9 for QPSK)
    dci.RedundancyVersion = 0;                          % Redundancy version (0,1,2,3)
    dci.SystemInformationIndicator = 0;                 % System information indicator. 0 for SIB1 (TS 38.212 Table 7.3.1.2.1-2)
    dci.ReservedBits = 0;
    
    % Use the same DCI for all PDCCH instances
    dci = repmat(dci,1,length(ss));

    % Configure as many control channels (PDCCH) as SS blocks in burst
    enabledPDCCH = transmittedBlocks & config.EnableSIB1;
    pdcch = getWavegenPDCCHConfig(ss,dci,enabledPDCCH);
    
    % Configure as many data channels (PDSCH) as SS blocks in burst
    ssSlotPeriod = vertcat(ss(:).SlotPeriodAndOffset);
    enablePDSCH = transmittedBlocks & config.EnableSIB1;
    pdsch = getWavegenPDSCHConfig(csetNRB,csetPattern,ssBurst.DMRSTypeAPosition,dci,ssSlotPeriod,enablePDSCH);
    
    % Configure waveform duration, frequency range, channel bandwidth
    numSubframes = 20;
    if any(scsSSB == [15 30])
        FR = 'FR1';
        BW = 100; % Max bandwidth for FR1
    else
        FR = 'FR2';
        BW = 400; % Max bandwidth for FR2
        if csetPattern == 1
            numSubframes = 25; % Increase waveform duration to ensure all 64 SSB + SIB1 fit
        end
    end
    
    % Map configuration to nrDLCarrierConfig
    wavegenConfig.SCSCarriers = carriers;
    wavegenConfig.BandwidthParts = bwps;
    wavegenConfig.SSBurst = ssBurst;
    wavegenConfig.CORESET = {coreset};
    wavegenConfig.SearchSpaces = mat2cell(ss(:),ones(1,length(ss)))';
    wavegenConfig.PDCCH = mat2cell(pdcch(:),ones(1,length(pdcch)))';
    wavegenConfig.PDSCH = mat2cell(pdsch(:),ones(1,length(pdsch)))';
    wavegenConfig.FrequencyRange = FR;
    wavegenConfig.ChannelBandwidth = BW;
    wavegenConfig.NumSubframes = numSubframes;
end

function [carrier,bwp] = getCommonSCSCarrierAndBWPConfig(ssb,minChanBW)
    
    scsCommon = ssb.SubcarrierSpacingCommon;
    
    carrier = nrSCSCarrierConfig;
    carrier.SubcarrierSpacing = scsCommon;
    carrier.NStartGrid = 0; 
    
    % Set the SCS carrier bandwidth to include CORESET 0
    msbIdx = floor(ssb.PDCCHConfigSIB1/16);
    scsSSB = hSSBurstSubcarrierSpacing(ssb.BlockPattern);
    scs = [scsSSB scsCommon];
    kssb = ssb.KSSB;
    [csetNRB,~,csetFreqOffset,~,csetTable] = hCORESET0Resources(msbIdx,scs,minChanBW,kssb);
    
    if isnan(csetNRB)
        idx = max(csetTable(~isnan(csetTable(:,2)),1));
        error('hSIB1WaveformConfiguration: For this SSB block pattern, common subcarrier spacing and minimum channel bandwidth, PDCCHConfigSIB1 must be lower than %d.', idx*16-1)
    end
    
    % CORESET 0 frequency offset from carrier center and carrier size
    c0 = (csetFreqOffset+10*scsSSB/scsCommon); 
    halfCarrierSize = max(c0,csetNRB-c0);

    % Adjust carrier bandwidth to account for CORESET being located in a
    % multiple of 6 RB from CRB0
    coff = mod(6-(halfCarrierSize+carrier.NStartGrid-c0),6);
    carrier.NSizeGrid = 2*(halfCarrierSize + coff);
    
    % Configure BWP for control and data. The carrier bandwidth above
    % ensures the BWP offset is defined in TS 38.213 Tables 13-1 to 13-10
    % and a multiple of 6 RB
    bwp = nrWavegenBWPConfig('BandwidthPartID',1,'Label','COMMON');
    bwp.SubcarrierSpacing = carrier.SubcarrierSpacing; 
    bwp.NSizeBWP = csetNRB; 
    bwp.NStartBWP = carrier.NStartGrid+carrier.NSizeGrid/2-c0;
    
end

function ssBurst = getWavegenSSBurstConfig(ssb)

    ssBurst = nrWavegenSSBurstConfig;
    ssBurst.BlockPattern = ssb.BlockPattern;
    ssBurst.TransmittedBlocks = ssb.TransmittedBlocks;
    
    % Update MIB content based on SIB1 config
    ssBurst.SubcarrierSpacingCommon = ssb.SubcarrierSpacingCommon;
    ssBurst.DMRSTypeAPosition = 3;

    % Select PDCCHConfigSIB1 to configure aspects related to the allocation
    % of CORESET 0 and PDCCH. The 4 MSB configure CORESET 0 resources
    % according to TS 38.213 Tables 13-1 through 13-10. The 4 LSB configure
    % Type0 PDCCH common search space monitoring occasions according to
    % Tables 13-11 through 13-15. The value of PDCCHConfigSIB1 below
    % configures CORESET pattern 1 with a single SS per slot
    ssBurst.PDCCHConfigSIB1 = 0*16+4; % (4MSB)*16 + (4LSB)     
    
    % Period of SS burst
    if ssBurst.SubcarrierSpacingCommon < 60 % FR1
        ssBurst.Period = 20; % ms
    else % FR2
        % Increase repetition period of the SSB for FR2 to ensure all
        % channels fit within one period
        ssBurst.Period = 40; % ms
    end
    
end

function SS = getType0CommonSearchSpaces(ssb,csetNRB,csetDuration,csetPattern)

    ssbTransmittedBlocks = ssb.TransmittedBlocks;
    L_max = length(ssbTransmittedBlocks);
    scs = [hSSBurstSubcarrierSpacing(ssb.BlockPattern) ssb.SubcarrierSpacingCommon];
    lsbIdx = mod(ssb.PDCCHConfigSIB1,16);
    maxAL = floor(log2(csetNRB*csetDuration/6))+1;
    slotsPerFrame = 10*scs(2)/15;
    
    framePeriod = ssb.Period/10; % 10 ms frames
    
    SS(1:L_max) = nrSearchSpaceConfig();
    ssID = 1;
    for ssbIdx = 0:L_max-1
        [ssSlot,ssFirstSym,~,frameOffset] = hPDCCH0MonitoringOccasions(lsbIdx,ssbIdx,scs,csetPattern,csetDuration);
        
        % Configure search space
        ss = nrSearchSpaceConfig();
        ss.CORESETID = 0;
        ss.SearchSpaceType = 'Common';
        ss.StartSymbolWithinSlot = ssFirstSym;
        
        period = framePeriod*slotsPerFrame; % Period in slots
        offset = mod(ssSlot + frameOffset*slotsPerFrame,period);
        ss.SlotPeriodAndOffset = [period offset];
        
        if csetPattern==1
            ss.Duration = 2;
        else % patterns 2 and 3
            ss.Duration = 1;
        end
        ss.NumCandidates = [0 0 4 2 1]; % TS 38.213 Table 10.1-1
        ss.NumCandidates(maxAL+1:end) = 0; % Limit aggregation levels to CORESET resources
        ss.SearchSpaceID = ssID;
        SS(ssID) = ss;
        ssID = ssID+1;
    end
end

% Configure SIB1 PDCCH. The length of the input search space ss, structure
% dci and enabled must be equal to the number of SSB.
function pdcch = getWavegenPDCCHConfig(ss,dci,enable)

    if length(dci) ~= length(ss)
        error('hSIB1WaveformConfiguration: The length of the DCI structure array must be equal to that of the search space configuration array.');
    end
    
    pdcch(1:length(ss)) = nrWavegenPDCCHConfig;

    for ch = 1:length(dci)
        dcibits = toBits(dci(ch));

        pdcch(ch).SearchSpaceID = ss(ch).SearchSpaceID;
        pdcch(ch).RNTI = 65535; % SI-RNTI. nRNTI = 0 as per TS 38.211 Section 7.3.2.3
        pdcch(ch).DataSource = dcibits;
        pdcch(ch).DataBlockSize = dci.Width;
        pdcch(ch).DMRSScramblingID = []; % Use NCellID TS 38.211 Section 7.3.2.3

        % Allocate PDCCH0 using Search Space instance
        pdcch(ch).Period = ss(ch).SlotPeriodAndOffset(1);
        pdcch(ch).SlotAllocation = ss(ch).SlotPeriodAndOffset(2);

        pdcch(ch).AggregationLevel = 2^(find(ss(1).NumCandidates,1,'last')-1);
        pdcch(ch).AllocatedCandidate = 1;

        % Enable PDCCH instances
        pdcch(ch).Enable = enable(ch);
    end

end

% Configure SIB1 PDSCH. The length of the input dci, slotPeriod and enabled
% must be equal to the number of SSB.
function pdsch = getWavegenPDSCHConfig(csetNRB,csetPattern,dmrsTypeAPosition,dci,ssSlotPeriod,enable)
    
    nch = length(dci);
    pdsch(1:nch) = nrWavegenPDSCHConfig;
    for ch = 1:nch

        [pdschTmp,K_0] = hSIB1PDSCHConfiguration(dci(ch),csetNRB,dmrsTypeAPosition,csetPattern);
        
        % Copy common PDSCH configuration to wavegen object
        pdsch(ch).RNTI = pdschTmp.RNTI;
        pdsch(ch).ReservedPRB = pdschTmp.ReservedPRB;
        pdsch(ch).Modulation = pdschTmp.Modulation;
        pdsch(ch).NumLayers = pdschTmp.NumLayers;
        pdsch(ch).MappingType = pdschTmp.MappingType;
        pdsch(ch).SymbolAllocation = pdschTmp.SymbolAllocation;
        pdsch(ch).PRBSet = pdschTmp.PRBSet;
        pdsch(ch).VRBToPRBInterleaving = pdschTmp.VRBToPRBInterleaving;
        pdsch(ch).NID = pdschTmp.NID;
        pdsch(ch).RNTI = pdschTmp.RNTI;
        pdsch(ch).DMRS = pdschTmp.DMRS;
        pdsch(ch).EnablePTRS = pdschTmp.EnablePTRS;

        % PDSCH configuration specific to the waveform generator includes
        % slot-wise allocation and transport layer aspects
        pdsch(ch).Period = ssSlotPeriod(ch,1);
        pdsch(ch).SlotAllocation = ssSlotPeriod(ch,2)+K_0;
        pdsch(ch).TargetCodeRate = hMCS(dci(ch).ModulationCoding);
        pdsch(ch).XOverhead = 0; % TS 38.214 Section 5.1.3.2
        pdsch(ch).RVSequence = dci(ch).RedundancyVersion;
        pdsch(ch).Enable = enable(ch);
        
    end
    
end

function RIV = hRIV(NSizeBWP,RBStart,LRBS)

    if LRBS<1 || LRBS>(NSizeBWP-RBStart)
        error('The number PDSCH allocation is out of the BWP limits.');
    end
    
    if (LRBS-1) <= floor(NSizeBWP/2)
        RIV = NSizeBWP*(LRBS-1)+RBStart;
    else
        RIV = NSizeBWP*(NSizeBWP-LRBS+1) + (NSizeBWP-1-RBStart);
    end
    
end