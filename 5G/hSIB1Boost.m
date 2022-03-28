%hSIB1Boost Boost SNR of SSB and associated SIB1 slot

%   Copyright 2020 The MathWorks, Inc.

function wave = hSIB1Boost(wave,waveConfig,waveInfo,ssbIdx,dBBoost)
    
    ssb = waveConfig.SSBurst;
    numBlocks = length(ssb.TransmittedBlocks);
    if ~(ssbIdx < numBlocks)
        error('hSIB1Boost: The specified SS block index (%d) must be lower than %d for a %s block pattern.',ssbIdx,numBlocks,ssb.BlockPattern);
    end
    
    if ~any( ssbIdx == (find(ssb.TransmittedBlocks)-1) )
        warning('hSIB1Boost: SS block specified is not active in burst.')
    end
    
    carriers = [waveConfig.SCSCarriers{:}];
    scs = [carriers.SubcarrierSpacing];
    ssbSCS = hSSBurstSubcarrierSpacing(ssb.BlockPattern);
    
    commonOfdmInfo = waveInfo.ResourceGrids.Info;
    ssbOfdmInfo = nrOFDMInfo(carriers(scs==ssbSCS).NSizeGrid,ssbSCS,'SampleRate',commonOfdmInfo.SampleRate);
    
    % Calculate the samples of the SSB to boost
    symLen = ssbOfdmInfo.SymbolLengths;
    symPerSub = length(symLen);
    samplesPerSub = ssbOfdmInfo.SampleRate*1e-3;
    
    burstStartSymbols = hSSBurstStartSymbols(ssb.BlockPattern,numBlocks);
    ssbSym = burstStartSymbols(ssbIdx+1);
    subframeOffset = floor(ssbSym/symPerSub)*samplesPerSub;
    symbolOffset = sum(symLen(1:mod(ssbSym,symPerSub)));
    ssblen = sum(symLen(1:mod(ssbSym+4,symPerSub))) - symbolOffset + 1;
    ssbFirstSample = subframeOffset + symbolOffset + 1;
    ssbFirstSample = ssbFirstSample + 0:ssb.Period*samplesPerSub:length(wave);
    boostedSamplesSSB = ssbFirstSample + (0:ssblen-1)';
    
    % Boost SNR of SSB symbols
    wave(boostedSamplesSSB) = db2mag(dBBoost)*wave(boostedSamplesSSB);
    
    if waveConfig.PDCCH{ssbIdx+1}.Enable
        % Boost SNR of associated SIB1 PDCCH and PDSCH slot
        % Calculate the samples of the SIB1 slot to boost
        symPerSlot = commonOfdmInfo.SymbolsPerSlot;
        slotPerSub = commonOfdmInfo.SlotsPerSubframe;
        samplesPerSub = commonOfdmInfo.SampleRate*1e-3;
        
        symLen = commonOfdmInfo.SymbolLengths;
        symPerSub = length(symLen);
        
        slot = waveConfig.PDCCH{ssbIdx+1}.SlotAllocation(1);
        subframeOffset = floor(slot/slotPerSub)*samplesPerSub;
        slotOffset = sum(symLen(1:mod(slot*symPerSlot,symPerSub)));
        slotFirstSample = subframeOffset + slotOffset + 1;
        slotFirstSample = slotFirstSample + 0:ssb.Period*samplesPerSub:length(wave);
        boostedSamplesSIB  = slotFirstSample + (0:sum(symLen(1:symPerSlot))-1)';
        
        % Boost SIB1 slot only if SIB1 samples have not been boosted yet.
        % For CORESET patterns 2 and 3, the SIB1 slot matches the location
        % of the SSB
        x = boostedSamplesSSB([1,end]);
        y = boostedSamplesSIB([1,end]);
        if ~(x(2) >= y(1) && y(2) >= x(1)) % No intersection of boosted samples
            wave(boostedSamplesSIB) = db2mag(dBBoost)*wave(boostedSamplesSIB);
        end
    end
end