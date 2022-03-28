% hSSBurstSubcarrierSpacing subcarrier spacing of SS block pattern

%   Copyright 2020 The MathWorks, Inc.

function scs = hSSBurstSubcarrierSpacing(blockPattern)

    if (strcmpi(blockPattern,'Case A'))
        scs = 15;
    elseif (any(strcmpi(blockPattern,{'Case B','Case C'})))
        scs = 30;
    elseif (strcmpi(blockPattern,'Case D'))
        scs = 120;
    elseif (strcmpi(blockPattern,'Case E'))
        scs = 240;
    end

end