% BitField Bit Field
%   F = BitField(WS,W,V) creates a bit field object F, with optional set of
%   expected bit widths WS, optional initial bit width W, and optional 
%   initial value V. By default, WS is [], indicating no restriction on the
%   field width, and W is 1 bit or min(WS) if WS is not empty. The field 
%   value defaults to 0 if not provided.
% 
%   BitField properties:
% 
%   Value       - Integer field value (value can be assigned from a bit vector)
%   Width       - Bit width of field
%   WidthValues - Set of expected, nominal bit widths for the field (Read-only)
%
%   BitField methods:
%
%   BitField - Class constructor
%   toBits   - Map the field value to information bit vector
%   fromBits - Map the field value from information bit vector
%   info     - Get information about the bit field
%   
%   Example:
%   % Create a BitField object where the field bit width is expected to be
%   % 0, 1, 5, or 10 bits. The initial width will be 0 (smallest of the set).
%   bf = BitField([0 1 5 10]);
% 
%   % Change width after construction to 5 bits.
%   bf.Width = 5;
% 
%   % Set the field value = 15, using a bit vector and then an integer.
%   bf.Value = [0 1 1 1 1];
%   bf.Value = 15;
%   bf.Value = hex2dec('F');
%   
%   % Map field value to bit vector.
%   bits = toBits(bf);
% 
%   See also MessageFormat.

%   Copyright 2021 The MathWorks, Inc.

classdef BitField
    
    properties (Dependent)
        Value   % Integer value of field
        Width   % Bit width of field
    end
    
    properties (SetAccess = private)
        WidthValues = [];   % Set of nominal widths defined for this field (Read-only)
    end
   
    methods
        
        % Constructor
        function obj = BitField(widthvalues,initwidth,initvalue)
        %BitField Construct an instance of BitField
        %   F = BitField(WS,W,V) creates a bit field object F, with optional set of
        %   expected bit widths WS, optional initial bit width W, and optional 
        %   initial value V. By default, WS is [], indicating no restriction on the
        %   field width, and W is 1 bit or min(WS) if WS is not empty. The field 
        %   value defaults to 0 if not provided.
            if nargin > 0
                obj.WidthValues = widthvalues;
                if nargin < 2
                    initwidth = min(widthvalues);  % If initial width is not supplying then use smallest nominal width
                end
                obj.Width = initwidth;
                if nargin > 2
                    obj.Value = initvalue;
                end       
            end          
        end
        
        function ifo = info(obj,opts)
        % info Get information about the field (value, width, width set)
            if nargin > 1 && any(strcmpi(opts,{'width','fieldsizes'}))
                ifo = obj.Width;
            elseif nargin > 1 && any(strcmpi(opts,'widthvalues'))
                ifo = obj.WidthValues;           
            else
                ifo = obj.Value;
            end
        end
               
        function bits = toBits(obj)
        % toBits Map field value into information bit vector (MSB maps to first bit in vector)

            % The most significant bit of each field is mapped to the lowest order information bit (i.e. first) of that field
            bits = int8(dec2bin(obj.Value,obj.Width) == '1')';
        end
    
        function obj = fromBits(obj,B)
        % fromBits Turn information bit vector into field value (first bit in vector maps to MSB of value)
            wl = min(obj.Width,length(B));
            obj.Value = sum(reshape(B(1:wl)~=0,1,[]).*2.^(wl-1:-1:0));
        end   
           
        % Set the field value from an integer or bit vector (modulo value by the field width)
        function obj = set.Value(obj,B)
            if length(B) > 1
                wl = min(obj.Width,length(B));
                B = sum(reshape(B(1:wl)~=0,1,[]).*2.^(wl-1:-1:0));
            end
            
            if ~isempty(B) || obj.Width == 0
                obj.IValue = mod(B,2^obj.Width);
            end 
        end
        
        % Get the integer field value (empty if the bit width is 0)
        function v = get.Value(obj)
            v = obj.IValue;
            if obj.IWidth == 0
               v = []; 
            end
        end
        
        % Set the field width
        function obj = set.Width(obj,B)
            if isempty(B) || B < 0
               error("The field bit width must be positive.") 
            end
            if ~isempty(obj.WidthValues) && ~any(B == obj.WidthValues)
                warning("The field bit width (%d) is not one of the set specified (%s) for this field.", B, strjoin(string(obj.WidthValues),','));
            end
            obj.IWidth = B;
            obj.Value = obj.Value;      % Adjust stored value wrt the new bit width
            if isempty(obj.Value) && B > 0
                obj.Value = 0;
            end
        end
        
        % Get the field width
        function w = get.Width(obj)
            w = obj.IWidth;
        end
        
    end
    
    properties (Access = private)
        IValue = 0;
        IWidth = 1;   
    end
   
end



