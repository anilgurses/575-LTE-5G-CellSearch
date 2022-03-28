%PDCCH_ConfigSIB1 PDCCH configuration for SIB1 information element in MIB
%   M = PDCCH_ConfigSIB1 creates a radio resource control (RRC) information
%   element of the physical downlink control channel (PDCCH) configuration
%   for the first system information block (SIB1) transmitted in the master
%   information block (MIB). PDCCH_ConfigSIB1 corresponds to the
%   information element PDCCH-ConfigSIB1, as defined in TS 38.331. All
%   PDCCH_ConfigSIB1 fields map, in order of the format definition, onto a
%   set of information bits. The mapping is such that the most significant
%   bit of each field is mapped to the lowest-order information bit for
%   that field. The number of bits associated with each field is fixed. The
%   information bits are encoded by the broadcast channel (BCH) and carried
%   on the physical broadcast channel (PBCH) as part of the MIB.
%
%   PDCCH_ConfigSIB1 properties (configurable):
%
%   controlResourceSetZero  - CORESET 0 configuration index
%   searchSpaceZero         - Search space 0 configuration index
%
%   PDCCH_ConfigSIB1 properties (read-only):
%
%   Width                   - Bit width of message format
%
%   Example:
%   % Create a PDCCH_ConfigSIB1 RRC information element and map the field
%   % values to information bits.
%
%   m = PDCCH_ConfigSIB1;
%   m.controlResourceSetZero = 0;
%   m.searchSpaceZero = 4;
% 
%   disp(toBits(m))
%
% See also MIB, MessageFormat, BitField.

%   Copyright 2021 The MathWorks, Inc.

classdef PDCCH_ConfigSIB1 < MessageFormat
    
    properties
        
        %controlResourceSetZero CORESET#0 configuration index
        % 4 most significant bits of PDCCH-ConfigSIB1 in the master
        % information block (MIB). controlResourceSetZero specifies the row
        % index of TS 38.213 Tables 13-1 to 13-10.
        controlResourceSetZero = BitField(4);
        
        %searchSpaceZero SearchSpace#0 configuration index
        % 4 least significant bits of PDCCH-ConfigSIB1 in the master
        % information block (MIB). searchSpaceZero specifies the row index
        % of TS 38.213 Tables 13-11 through 13-15 for the PDCCH monitoring
        % occasions.
        searchSpaceZero = BitField(4);
        
    end
    
end
