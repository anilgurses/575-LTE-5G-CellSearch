%hPDSCHTimeAllocationTables PDSCH time domain resource allocation tables

%   Copyright 2020-2021 The MathWorks, Inc.

function out = hPDSCHTimeAllocationTables()

    % K_0: slot offset from DCI slot scheduling PDSCH
    % S: Start symbol of allocated PDSCH
    % L: Length in symbols of allocated PDSCH
    
    persistent restables;
    
    if (isempty(restables))
    
        dup = @(x)reshape(repmat(x,2,1),[],1);
        rep = @(x,n)repmat(x.',n,1);

        % TS 38.214 Table 5.1.2.1.1-2
        rowIndex = dup(1:16);
        DMRSTypeAPosition = rep([2 3],16);
        %                  row   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
        PDSCHMappingType = dup(["A" "A" "A" "A" "A" "B" "B" "B" "B" "B" "B" "A" "A" "A" "B" "B"]);
        K_0 = rep([0 0],16);
        % row  1  1  2  2  3  3  4  4  5  5  6  6  7  7          8  9 10 11 12 13 14 15 16
        S = [[ 2  3  2  3  2  3  2  3  2  3  9 10  4  6].'; dup([5  5  9 12  1  1  2  4  8])];
        L = [[12 11 10  9  9  8  7  6  5  4  4  4  4  4].'; dup([7  2  2  2 13  6  4  7  4])];
        defaultA = table(rowIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

        % TS 38.214 Table 5.1.2.1.1-4
        rowIndex = [dup(1:15); 16];
        DMRSTypeAPosition = [rep([2 3],15); NaN];
        %                   row   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15     16 
        PDSCHMappingType = [dup(["B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "A" "B"]); ""];
        %      row  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15    16
        K_0 = [dup([0  0  0  0  0  1  1  0  0  0  0  0  0  0  1]); NaN];
        %    row  1  2  3  4  5  6  7  8  9 10 11 12 13    14  14  15  15  16
        S = [dup([2  4  6  8 10  2  4  2  4  6  8 10  2]);  2;  3;  2;  2; NaN];
        L = [dup([2  2  2  2  2  2  2  4  4  4  4  4  7]); 12; 11;  4;  4; NaN];
        defaultB = table(rowIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

        % TS 38.214 Table 5.1.2.1.1-5
        rowIndex = [dup(1:5); 6; 7; dup(8:16)];
        DMRSTypeAPosition = [rep([2 3],5); NaN; NaN; rep([2 3],9)];
        %                   row   1   2   3   4   5     6   7         8   9  10  11  12  13  14  15  16
        PDSCHMappingType = [dup(["B" "B" "B" "B" "B"]); ""; ""; dup(["B" "B" "B" "B" "B" "B" "A" "A" "A"])];
        K_0 = [rep([0 0],5); NaN; NaN; rep([0 0],9)];
        %    row  1  2  3  4  5     6    7        8  9 10 11 12 13    14  14      15 16
        S = [dup([2  4  6  8 10]); NaN; NaN; dup([2  4  6  8 10  2]);  2;  3; dup([0  2])];
        L = [dup([2  2  2  2  2]); NaN; NaN; dup([4  4  4  4  4  7]); 12; 11; dup([6  6])];
        defaultC = table(rowIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

        restables = {defaultA defaultB defaultC};
        
    end
    
    out = restables;
    
end