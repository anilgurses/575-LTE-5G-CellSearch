%hPDCCH0Configuration PDCCH configuration in Type0 CORESET/SS
%   PDCCH = hPDCCH0Configuration(SSBINDEX,INITSYSTEMINFO,SCS,NCELLID,MINCHANBW)
%   returns a Type 0 physical downlink control channel configuration object
%   PDCCH configured by synchronization signal block (SSB) index SSBINDEX,
%   initial system information INITSYSTEMINFO, SSB and common subcarrier
%   spacing SCS ([SSB,COMMON]), cell identity NCELLID and minimum channel
%   bandwidth MINCHABW.

%   Copyright 2020-2021 The MathWorks, Inc.

function pdcch = hPDCCH0Configuration(ssbIndex,initSystemInfo,scs,ncellid,minChanBW)

    % Obtain CORESET configuration parameters based on PDCCH-ConfigSIB1 in
    % MIB. See TS 38.213 Section 13
    msbIdx = initSystemInfo.PDCCHConfigSIB1.controlResourceSetZero;
    kssb = initSystemInfo.k_SSB;
    [csetNRB,csetDuration,~,csetPattern] = hCORESET0Resources(msbIdx,scs,minChanBW,kssb);
    
    coreset = nrCORESETConfig();
    coreset.FrequencyResources = ones(1,csetNRB/6);
    coreset.Duration = csetDuration;
    coreset.CCEREGMapping = 'interleaved'; % 38.211-7.3.2.2
    coreset.REGBundleSize = 6; 
    coreset.InterleaverSize = 2; 
    coreset.ShiftIndex = ncellid;
    coreset.CORESETID = 0;
    
    % Get search space monitoring parameters based on PDCCH-ConfigSIB1 in MIB
    % TS 38.213 Section 13
    lsbIdx = initSystemInfo.PDCCHConfigSIB1.searchSpaceZero;
    [ssSlot,ssFirstSym,isOccasion] = hPDCCH0MonitoringOccasions(lsbIdx,ssbIndex,scs,csetPattern,csetDuration,initSystemInfo.NFrame);
    
    % If there are no PDCCH monitoring occasions in this SFN, check the next
    slotsPerFrame = 10*initSystemInfo.SubcarrierSpacingCommon/15;
    if ~isOccasion
        [ssSlot,ssFirstSym,isOccasion] = hPDCCH0MonitoringOccasions(lsbIdx,ssbIndex,scs,csetPattern,csetDuration,initSystemInfo.NFrame+1);
        ssSlot = ssSlot+slotsPerFrame;
    end

    assert(isOccasion,'This frame does not contain any PDCCH monitoring occasion.');

    % Configure search space
    ss = nrSearchSpaceConfig();
    ss.CORESETID = coreset.CORESETID;
    ss.SearchSpaceType = 'Common';
    ss.StartSymbolWithinSlot = ssFirstSym;
    if csetPattern==1
        ss.SlotPeriodAndOffset = [2*slotsPerFrame ssSlot];
        ss.Duration = 2;
    else % patterns 2 and 3
        ss.SlotPeriodAndOffset = [2*slotsPerFrame ssSlot];
        ss.Duration = 1;
    end
    ss.NumCandidates = [0 0 4 2 1]; % TS 38.213 Table 10.1-1
    
    % Limit aggregation levels to CORESET resources
    maxAL = floor(log2(csetNRB*csetDuration/6))+1;
    ss.NumCandidates(maxAL+1:end) = 0; 
    
    % According to TS 38.212 Section 7.3.1.0, for CORESET 0 the BWP is the
    % size of the CORESET.
    pdcch = nrPDCCHConfig();
    pdcch.NStartBWP = 0;
    pdcch.NSizeBWP = csetNRB;
    pdcch.CORESET = coreset;
    pdcch.SearchSpace = ss;
    pdcch.RNTI = 0; % TS 38.211 Section 7.3.2.3
    pdcch.DMRSScramblingID = ncellid; % TS 38.211 Section 7.3.2.3
    
end