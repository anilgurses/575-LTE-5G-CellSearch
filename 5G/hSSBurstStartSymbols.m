%hSSBurstStartSymbols Starting OFDM symbols of SS burst

%   Copyright 2020 The MathWorks, Inc.

function ssbStartSymbols = hSSBurstStartSymbols(ssbBlockPattern,Lmax)

    % 'alln' gives the overall set of SS block indices 'n' described in 
    % TS 38.213 Section 4.1, from which a subset is used for each Case A-E
    alln = [0; 1; 2; 3; 5; 6; 7; 8; 10; 11; 12; 13; 15; 16; 17; 18];
    
    cases = {'Case A' 'Case B' 'Case C' 'Case D' 'Case E'};
    m = [14 28 14 28 56];
    i = {[2 8] [4 8 16 20] [2 8] [4 8 16 20] [8 12 16 20 32 36 40 44]};
    nn = [2 1 2 16 8];
    
    caseIdx = find(strcmpi(ssbBlockPattern,cases));
    if (any(caseIdx==[1 2 3]))
        if (Lmax==4)
            nn = nn(caseIdx);
        elseif (Lmax==8)
            nn = nn(caseIdx) * 2;
        else
            error('For %s, the SSBTransmitted bitmap must be of length 4 or 8.',cases{caseIdx});
        end
    else
        if (Lmax==64)
            nn = nn(caseIdx);
        else
            error('For %s, the SSBTransmitted bitmap must be of length 64.',cases{caseIdx});
        end
    end
    
    n = alln(1:nn);
    ssbStartSymbols = (i{caseIdx} + m(caseIdx)*n).';
    ssbStartSymbols = ssbStartSymbols(:).';
    
end