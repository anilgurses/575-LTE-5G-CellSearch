%hPDCCH0MonitoringOccasions Type 0 PDCCH monitoring occasions in CSS
%   [N0,FSYM,ISOCASSION,FRAMEOFF] = ...
%   hPDCCH0MonitoringOccasions(INDEX,SSBIDX,SCS,CSETPAT,CSETDURATION,SFN)
%   returns the slot index N0 and first symbol index FSYM from TS 38.213
%   Section 13 Tables 13-11 through 13-15 for the PDCCH monitoring
%   occasions corresponding to a row index INDEX (4 LSB of PDCCH-ConfigSIB1
%   in MIB). The table is selected using the two-element vector SCS
%   ([SSB,COMMON]), coreset pattern CSETPAT and coreset duration
%   CSETDURATION. The output ISOCCASION indicates if there are control
%   monitoring occasions associated to the SSB index SSBIDX in the system
%   frame number specified in SFN. The output FRAMEOFF indicates the frame
%   offset from SFN of the monitoring occasion associated to N0.

%   Copyright 2020-2021 The MathWorks, Inc.

function [n0,fSym,isOccasion,frameOff] = hPDCCH0MonitoringOccasions(index,ssbIdx,scs,pat,duration,sfn)

    scsSSB = scs(1);
    scsCommon = scs(2);
    symPerSlot = 14;
    
    if pat == 1      
        if (scsSSB==15 || scsSSB==30) % Frequency Range 1
            if mod(ssbIdx,2) == 0
                fSym = 0;
            else
                fSym = duration;
            end
            
            tab = [ 0,    1, 2,    3, 4,    5, 6,    7, 8, 9, 10, 11, 12, 13, 14, 15;...
                    0,    0, 2,    2, 5,    5, 7,    7, 0, 5,  0,  0,  2,  2,  5,  5;...
                    1,    2, 1,    2, 1,    2, 1,    2, 1, 1,  1,  1,  1,  1,  1,  1;...
                    1,  1/2, 1,  1/2, 1,  1/2, 1,  1/2, 2, 2,  1,  1,  1,  1,  1,  1;...
                    0, fSym, 0, fSym, 0, fSym, 0, fSym, 0, 0,  1,  2,  1,  2,  1,  2];
        else % Frequency Range 2
            if mod(ssbIdx,2) == 0 % Even SSB index
                fSym1 = 0;
                fSym2 = 0;
            else
                fSym1 = 7;
                fSym2 = duration;
            end
            
            tab = [ 0,     1,   2,     3, 4,     5,     6,     7,     8,    9,    10,    11, 12, 13, 14, 15;...
                    0,     0, 2.5,   2.5, 5,     5,     0,   2.5,     5,  7.5,   7.5,   7.5,  0,  5, NaN(1,2);...
                    1,     2,   1,     2, 1,     2,     2,     2,     2,    1,     2,     2,  1,  1, NaN(1,2);...
                    1,   1/2,   1,   1/2, 1,   1/2,   1/2,   1/2,   1/2,    1,   1/2,   1/2,  2,  2, NaN(1,2);...
                    0, fSym1,   0, fSym1, 0, fSym1, fSym2, fSym2, fSym2,    0, fSym1, fSym2,  0,  0, NaN(1,2)];
        end
        
        % Extract info from table
        t = tab(:,index==tab(1,:));
        O = t(2);
        M = t(4);
        fSym = t(5);
        
        mu = log2(scsCommon/15);
        slotsPerFrame = 10*2^mu;
        
        % Slot of monitoring occasion
        slot = O*2^mu + floor(ssbIdx*M);
        frameOff = floor(slot/slotsPerFrame);
        n0 = mod(slot,slotsPerFrame);
        
        isOccasion = true;
        if nargin == 6
            if ~isnan(frameOff) && ~xor(mod(frameOff,2),mod(sfn,2))
                isOccasion = true;
            else
                isOccasion = false;
            end
        end
        
    else
        if pat == 2
            if all(scs == [120 60])
                ssbStartSymbols = hSSBurstStartSymbols('Case D',64);
                nc = floor(ssbStartSymbols(ssbIdx+1)/(symPerSlot*scsSSB/scsCommon));
                fSymTable = [0 1 6 7];
            elseif all(scs == [240 120])
                ssbStartSymbols = hSSBurstStartSymbols('Case E',64);
                nc = floor(ssbStartSymbols(ssbIdx+1)/(symPerSlot*scsSSB/scsCommon));
                if any(mod(ssbIdx,8)==[4 5])
                    nc = nc-1;
                end
                fSymTable = [0 1 2 3 12 13 0 1];
            else
                error('For CORESET pattern 2, common subcarrier spacing must be 60 or 120 kHz for SS block pattern ''Case D'' or ''Case E'', respectively.')
            end
            fSym = fSymTable(mod(ssbIdx,length(fSymTable))+1);
            
        elseif pat == 3
            if all(scs == [120 120])
                ssbStartSymbols = hSSBurstStartSymbols('Case D',64);
                nc = floor(ssbStartSymbols(ssbIdx+1)/(symPerSlot*scsSSB/scsCommon));
                fSymTable = [4 8 2 6];
                fSym = fSymTable(mod(ssbIdx,length(fSymTable))+1);
            else
                error('For CORESET pattern 3, common subcarrier spacing must be 120 kHz and SS block pattern must be ''Case D''.')
            end
        else
            error('CORESET pattern not supported.')
        end
        isOccasion = true;
        frameOff = 0;
        n0 = nc;
        
    end
    
end