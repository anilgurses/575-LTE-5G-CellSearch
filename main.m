addpath('./LTE');
addpath('./5G');


mode = "LTE";

if mode == "LTE"
    lteCellSearchSync(17);
else
    NRCellSearchSync(102);
end
