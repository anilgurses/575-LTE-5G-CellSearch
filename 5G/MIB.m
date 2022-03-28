%MIB MIB RRC message transmitted on BCH
%   M = MIB creates a master information block (MIB) radio resource control
%   (RRC) message that contains the MIB bit fields, as  defined in TS
%   38.331. The message enables the encoding and decoding of the MIB RRC
%   message bits from and to message field values. All the fields map, in
%   order of the format definition, onto a set of information bits. The
%   mapping is such that the most significant bit of each field is mapped
%   to the lowest-order information bit for that field. The number of bits
%   associated with each field is fixed. The information bits are encoded
%   by the broadcast channel (BCH) and carried on the physical broadcast
%   channel (PBCH). The BCH transport block is the RRC message
%   BCCH-BCH-Message, consisting of a leading 0 bit and 23 bits
%   corresponding to the MIB. The leading bit signals the message type
%   transmitted (MIB or empty sequence).
%
%   MIB properties (configurable):
%
%   systemFrameNumber       - System frame number
%   subCarrierSpacingCommon - Subcarrier spacing for SIB1 in PDSCH
%   ssb_SubcarrierOffset    - Subcarrier offset kSSB
%   dmrs_TypeA_Position     - First DM-RS symbol position for SIB1 in PDSCH
%   pdcch_ConfigSIB1        - PDCCH-ConfigSIB1 information element for downlink control configuration
%   cellBarred              - Cell barred
%   intraFreqReselection    - Intra frequency reselection
%   spare                   - Spare bit
%
%   MIB properties (read-only):
%
%   Width                   - Bit width of message format
%
%   Example:
%   % Create a MIB RRC message and map the field values to information bits.
%
%   m = MIB;
%   m.systemFrameNumber = 1;
%   m.subCarrierSpacingCommon = 0;                  % SCS 15 for FR1
%   m.ssb_SubcarrierOffset = 4;                     % kSSB = 4 for FR1
%   m.dmrs_TypeA_Position = 0;                      % DM-RS position 2
%   m.pdcch_ConfigSIB1.controlResourceSetZero = 0;
%   m.pdcch_ConfigSIB1.searchSpaceZero = 4;
%   m.cellBarred = 0;                               % Cell barred
%   m.intraFreqReselection = 0;                     % Reselection allowed
% 
%   disp(toBits(m))
%
% See also PDCCH_ConfigSIB1, MessageFormat, BitField.

%   Copyright 2021 The MathWorks, Inc.

classdef MIB < MessageFormat
    
    properties

        %systemFrameNumber System frame number 
        % 6 most significant bits (MSB) of the 10-bit system frame number
        % (SFN). The 4 LSB of the SFN are conveyed in the BCH transport
        % block as part of channel coding (outside the MIB encoding).
        systemFrameNumber = BitField(6); 

        %subCarrierSpacingCommon Subcarrier spacing for SIB1 in PDSCH
        % The value 0 (scs15or60) corresponds to 15 kHz and the value 1
        % (scs30or120) corresponds to 30 kHz when the UE acquires this MIB
        % on an frequency range 1 (FR1) carrier frequency. The value 0
        % (scs15or60) corresponds to 60 kHz and the value 1 (scs30or120)
        % corresponds to 120 kHz when the UE acquires this MIB on an FR2
        % carrier frequency.
        subCarrierSpacingCommon = BitField(1);
        
        %ssb_SubcarrierOffset Subcarrier offset kSSB
        % Frequency domain offset between SSB and the overall resource
        % block grid in number of subcarriers. The value range of this
        % field may be extended by an additional most significant bit
        % encoded within PBCH as specified in TS 38.213.
        ssb_SubcarrierOffset = BitField(4);

        %dmrs_TypeA_Position First DM-RS symbol position for SIB1 in PDSCH
        % PDSCH DM-RS symbol position bit. The value 0 indicates DM-RS
        % symbol position 2 and the value 1 indicates position 3.
        dmrs_TypeA_Position = BitField(1);

        %pdcch_ConfigSIB1 PDCCH-ConfigSIB1 information element for downlink control configuration
        % RRC information element containing controlResourceSetZero and
        % searchSpaceZero configuration.
        pdcch_ConfigSIB1 = PDCCH_ConfigSIB1;

        %cellBarred Cell barred
        % cellBarred configures permission for UEs to camp in the cell. A
        % false value indicates that cell does not allow UEs to camp.
        cellBarred = BitField(1);
        
        %intraFreqReselection Intra frequency reselection
        % intraFreqReselection configures cell selection or reselection to
        % intra-frequency cells when the highest ranked cell is barred, or
        % treated as barred by the UE. A false value indicates frequency
        % reselection is allowed.
        intraFreqReselection = BitField(1);

        %spare Spare bit
        spare = BitField(1);
        
    end

end
