%DCIFormat1_0_SIRNTI DCI format 1_0 with CRC scrambled by SI-RNTI
%   M = DCIFormat1_0_SIRNTI(NDLRB) creates a downlink control information
%   (DCI) format 1_0 message with CRC scrambled by the system information
%   radio network temporary identifier (SI-RNTI) for a CORESET 0 of size
%   NDLRB. The message contains the information described in TS 38.212
%   Clause 7.3.1.2.1. The message enables the encoding and decoding of the
%   DCI message bits from and to message field values. All the fields map,
%   in order of the format definition, onto a set of information bits. The
%   mapping is such that the most significant bit of each field is mapped
%   to the lowest-order information bit for that field. The information
%   bits are encoded and carried on the physical downlink control channel
%   (PDCCH).
%
%   M = DCIFormat1_0_SIRNTI(NDLRB,SHAREDSPECTRUM) enables to create a DCI
%   format 1_0 message with CRC scrambled by SI-RNTI with shared spectrum
%   enabled. A non-zero value of SHAREDSPECTRUM enables shared spectrum and
%   extends the number of reserved bits from 15 to 17.
%
%   DCIFormat1_0_SIRNTI properties (configurable):
%
%   FrequencyDomainResources    - Frequency domain resource assignment
%   TimeDomainResources         - Time domain resource assignment
%   VRBToPRBMapping             - VRB-to-PRB mapping
%   ModulationCoding            - Modulation and coding scheme
%   RedundancyVersion           - Redundancy version
%   SystemInformationIndicator  - System information indicator
%   ReservedBits                - Reserved bits
%
%   DCIFormat1_0_SIRNTI properties (read-only):
%
%   Width                       - Bit width of message format
%
%   Example:
%   % Create a DCIFormat1_0_SIRNTI message for a 48 RB CORESET 0 and map
%   % the field values to information bits.
%
%   NDLRB = 48;
%   dci = DCIFormat1_0_SIRNTI(NDLRB);
%   dci.FrequencyDomainResources = 1;
%   dci.ModulationCoding = 2;
%   disp(toBits(dci))
%
% See also MessageFormat, BitField.

%   Copyright 2021 The MathWorks, Inc.

classdef DCIFormat1_0_SIRNTI < MessageFormat
    
    properties

        %FrequencyDomainResources Frequency domain resource assignment
        % The bit width depends on the size of CORESET 0 provided to the
        % constructor, as described in TS 38.212 Clause 7.3.1.2.1.
        FrequencyDomainResources = BitField();

        %TimeDomainResources Time domain resource assignment
        % 4 bits as defined in TS 38.214 Subclause 5.1.2.1.
        TimeDomainResources = BitField(4);

        %VRBToPRBMapping VRB-to-PRB mapping
        % 1 bit according to TS 38.212 Table 7.3.1.2.2-5.
        VRBToPRBMapping = BitField(1);

        %ModulationCoding Modulation and coding scheme
        % 5 bits as defined in TS 38.214 Subclause 5.1.3.
        ModulationCoding = BitField(5);
        
        %RedundancyVersion Redundancy version
        % 2 bits as defined in TS 38.212 Table 7.3.1.1.1-2.
        RedundancyVersion = BitField(2);

        %SystemInformationIndicator System information indicator
        % 1 bit as defined in TS 38.212 Table 7.3.1.2.1-2.
        SystemInformationIndicator = BitField(1);

        %ReservedBits Reserved bits
        % 17 bits for operation in a cell with shared spectrum channel
        % access. Otherwise 15 bits.
        ReservedBits = BitField(15);

    end
    
    methods
       
        function obj = DCIFormat1_0_SIRNTI(nsizebwp,sharedspectrum)
            
            % NDLRB for CORESET 0
            N = ceil(log2(nsizebwp*(nsizebwp+1)/2));
            obj.FrequencyDomainResources = BitField(N);
            
            % If operation in a Release 16 cell with shared spectrum
            % channel access then adjust the field size
            if nargin>1 && sharedspectrum
                obj.ReservedBits = BitField(17);
            end

        end
    
    end
    
end
