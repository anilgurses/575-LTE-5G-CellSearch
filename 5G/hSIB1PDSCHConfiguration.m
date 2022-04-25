%hSIB1PDSCHConfiguration PDSCH configuration for SIB1
%   PDSCH = hSIB1PDSCHConfiguration(DCI,NSIZEBWP,DRMSTYPEAPOSITION,CSETPAT)
%   returns a physical downlink shared channel configuration object PDSCH
%   needed to demodulate a PDSCH carrying the system information block 1.
%   The PDSCH object is configured according to the downlink control
%   information structure DCI, size of the bandwidth part NSIZEBWP, PDSCH
%   DM-RS Type A position DMRSTYPEPOSITION and coreset pattern CSETPAT.

%   Copyright 2020-2021 The MathWorks, Inc.

function [pdsch,K_0] = hSIB1PDSCHConfiguration(dci,NSizeBWP,DMRSTypeAPosition,pat,modulation)

    pdsch = nrPDSCHConfig();
    pdsch.NSizeBWP = [];
    pdsch.NStartBWP = [];
    
    % Configure PDSCH from DCI message
    [Lrbs,RBstart] = hDecodeRIV(NSizeBWP,dci.FrequencyDomainResources);
    pdsch.PRBSet = RBstart + (0:(Lrbs-1));
    pdsch.VRBToPRBInterleaving = dci.VRBToPRBMapping;
    pdsch.RNTI = 65535; % SI-RNTI 
    
    % Select applicable PDSCH time domain resource allocation table based
    % on RNTI and CORESET pattern (TS 38.214 Table 5.1.2.1.1-1)
    % SI-RNTI Type0 Common
    restables = hPDSCHTimeAllocationTables(); 
    restable = restables{pat};
    restable = restable(restable.rowIndex==(dci.TimeDomainResources+1),:);
    resalloc = restable(restable.DMRSTypeAPosition==DMRSTypeAPosition,:);
    
    K_0 = resalloc.K_0;
    
    pdsch.MappingType = resalloc.PDSCHMappingType;
    pdsch.SymbolAllocation = [resalloc.S resalloc.L];
   
    % Configure PDSCH from MIB / SSB
    pdsch.NID = []; % Cell identity (TS 38.211 Section 7.3.1.1)
    pdsch.DMRS.DMRSTypeAPosition = DMRSTypeAPosition;
    pdsch.DMRS.NIDNSCID = []; % Cell identity (TS 38.211 Section 7.4.1.1)
    
    % Configure PDSCH from other relevant rules
    % TS 38.214 Section 5.1.3.1
    pdsch.Modulation = modulation;
    % TS 38.211 Section 7.4.1.1
    pdsch.DMRS.NSCID = 0;
    % TS 38.214 Section 5.1.6.2
    pdsch.NumLayers = 1; 
    pdsch.DMRS.DMRSPortSet = 0;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 1;
    L = pdsch.SymbolAllocation(2);
    if (L==2)
        pdsch.DMRS.NumCDMGroupsWithoutData = 1;
    else
        pdsch.DMRS.NumCDMGroupsWithoutData = 2;
    end
    if (strcmpi(pdsch.MappingType,'A'))
        pdsch.DMRS.DMRSAdditionalPosition = 2;
    else % 'B'
        switch L
            case {2,4}
                pdsch.DMRS.DMRSAdditionalPosition = 0;
            case 7
                pdsch.DMRS.DMRSAdditionalPosition = 1;
        end
    end
    
    pdsch.DMRS.DMRSReferencePoint = 'PRB0';
    
    pdsch.EnablePTRS = false; % TS 38.214 Section 5.1.6.3
    
end

function [Lrbs,RBstart] = hDecodeRIV(NSizeBWP,RIV)
    
    Lrbs = floor(RIV / NSizeBWP) + 1;
    
    RBstart = RIV - ((Lrbs - 1) * NSizeBWP);
    
    if (Lrbs > NSizeBWP - RBstart)
        
        Lrbs = NSizeBWP - Lrbs + 2;
        RBstart = NSizeBWP - 1 - RBstart;
        
    end
    
end
