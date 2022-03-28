% MessageFormat Message format base class
%   MessageFormat provides a base class to derive concrete message classes.
%   Use BitField properties in the derived classes to add format specific
%   message fields. Use the message objects to map fields into information
%   bit vectors, and from bit vectors back into message fields. The class
%   provides optional zero padding of the information bits for message
%   length alignment.
%   
%   BitField properties in derived classes store both field values and
%   bitwidths. MessageFormat overloads the subsref and subsasgn functions
%   to provide direct access to field value data. It also overloads the
%   display functionality to show the field values.
% 
%   MessageFormat properties:
% 
%   AlignedWidth - Target width of aligned message format (use [] for no automatic alignment padding)
%   Width        - Bit width of message format (including alignment padding) (Read-only)
%   PaddingWidth - Bit width of message padding (Read-only)
%        
%   MessageFormat methods:
%
%   toBits        - Map the message format to information bits
%   fromBits      - Map the message format from information bits
%   info          - Get information about message format
%   
%   Example:
%   % Define a message format consisting of a 1 bit field followed by a 
%   % 5 bit field.
%
%   classdef MyMessage < MessageFormat
%       properties
%           messageField1 = BitField(1);
%           messageField2 = BitField(5);
%       end
%     end
% 
%   See also BitField. 

%   Copyright 2021 The MathWorks, Inc.

classdef MessageFormat < matlab.mixin.CustomDisplay

    properties
        AlignedWidth = [];  % Target width of aligned message format (use [] for no automatic alignment padding)
    end

    properties (SetAccess = private, Dependent)    
        Width;              % Bit width of message format (including alignment padding) (Read-only)
        PaddingWidth;       % Bit width of message padding (Read-only)
    end

    properties (Dependent, Access = private)    
        ManagedPaddingField;  % Internal bit field that used to implement positive alignment padding
    end
        
    methods
        
        function bits = toBits(obj)
        % toBits Map message fields into information bit vector

            propName = getActiveFieldProperties(obj);  % Get field property names (inc padding) to use in mapping 
            bits = logical([]);
            for n = 1:numel(propName)
                bits = [bits; obj.(propName{n}).toBits];  %#ok<AGROW>
            end 
        end
        
        function [obj, bitpos] = fromBits(obj,bits)
        % fromBits Turn information bit vector into message field values

            % Check input length
            mwid = obj.Width;  % Format length
            if length(bits) < mwid
                bits(mwid) = 0;     % Resize/pad the input if required
            end

            [propName] = getActiveFieldProperties(obj); % Get field property names (inc padding) to use in mapping
            % Walk through the fields
            bitpos = 1;
            for n = 1:numel(propName)
                bitwid = obj.(propName{n}).Width;
                obj.(propName{n}) = obj.(propName{n}).fromBits(bits(bitpos:bitpos+bitwid-1));    % Call 'fromBits' on the individual fields
                bitpos = bitpos + bitwid;
            end

            bitpos = bitpos-1;  % Adjust, so that it reflects the number of bits read
        end
        
        function ifo = info(obj,opts)     
        % info Get information about the message
        
            if nargin < 2
                opts = [];
            end
            % Get all the relevant property names
            [propName,propNameOut] = getActiveFieldProperties(obj);
            
            % Create a structure to capture the info results
            propNameS = propNameOut';          % Create row of field names 
            propNameS{2,end} = [];             % Add second row in the cell array for the info 'values'
            ifo = struct(propNameS{:});        % Use this name/value list to create the info structure
            % Despatch the info options to properties and assignment the results
            for n = 1:numel(propName)          % Walk through all the fields
                ifo.(propNameOut{n}) = obj.(propName{n}).info(opts);
            end
        end
    
        % Calculate the alignment adjustment required
        function w = get.PaddingWidth(obj)
            if isempty(obj.AlignedWidth)
                w = 0;
            else
                w = obj.AlignedWidth - calcWidth(obj);  % This function call does not include any padding, so will not re-enter
            end
        end

        % Calculate the total bit width (including padding)
        function w = get.Width(obj)
            if  isempty(obj.AlignedWidth)
                w = calcWidth(obj);   % This function call does not include any padding
            else
                w = obj.AlignedWidth;
            end
        end
      
        % Get managed padding bit field
        function value = get.ManagedPaddingField(obj)
            value = BitField(max(0,obj.PaddingWidth));
        end

        % Set managed padding bit field (nop writeback)
        function obj = set.ManagedPaddingField(obj,~)
            % We can ignore any attempt to write back to the managed field (value = BitField(obj.PaddingWidth))
        end

        % Implement read-only semantics
        function obj = set.Width(obj,~) %#ok<INUSD>
            error("You cannot set the read-only property 'Width'");     
        end
        function obj = set.PaddingWidth(obj,~) %#ok<INUSD>
            error("You cannot set the read-only property 'PaddingWidth'");     
        end
        
        % Override of subscript assignment operation to allow direct assignment of BitField values 
        function obj = subsasgn(obj, S, B)
           
            % Examples of possible subscripting:
            % dci.DCIField
            % dci(x).DCIField
            % dci_.DCIField.BitFieldProps
            % dci_.MessageField.DCIField
            % dci_.MessageField.DCIField.BitFieldProps
            % 
            % Peek at type of value between subscripted, and if a BitField, extend subscripting 
            % to address the BitField's Value property, rather than the entire BitField object
            vpeek = builtin('subsref', obj, S);  
            if isa(vpeek,'BitField') && ~isa(B,'BitField')
                 S = [S,struct('type','.','subs','Value')];  % Extend subscripting with '.Value'
            end
            obj = builtin('subsasgn', obj, S, B);   % Pass on the subscript assignment operation  

        end
        
        % Override of subscript reference operation
        function V = subsref(obj, S)
                        
            % Staged lookup to trap any BitField properties
            V = builtin('subsref', obj, S);   % Pass on the subscript assignment operation 
            if isa(V,'BitField')
                V = V.Value;  
            end
            
        end
        
    end
    
    methods (Access = private)
        
        % Get field properties that will be treated a countable message fields 
        function [p,po] = getActiveFieldProperties(obj,excpadding)
            
            if nargin==1
                excpadding = false;  % If a second argument wasn't provided then include padding 
            end
         
            p = properties(obj);    % Start with all property names of the (derived) class
            p = p(1:end-3);         % These are the private facing property names (don't include the 3 public width related properties defined in this class)
            po = p;                 % These are the public facing property names

            % Include the managed padding field for size alignment if active?
            if ~(excpadding || isempty(obj.AlignedWidth))   % If padding requested to be included AND aligned width defined
               % Substitute or insert the managed padding field name into the property names
               paddingpos = find(strcmpi('padding',p),1); 
               if isempty(paddingpos)    % If no explicit 'padding' field found then append at end
                   paddingpos = length(p)+1;
               end
               p{paddingpos} = 'ManagedPaddingField';  % Internal facing field property name
               po{paddingpos} = 'Padding';             % External facing field property name
            else
                % Ensure no padding appears in list
               paddingpos = strcmpi('padding',p); 
               p(paddingpos) = [];     % Internal facing field property name
               po(paddingpos) = [];    % External facing field property name
            end
        end

    end    
    
    methods (Access = protected)
    
        function b = isInactiveProperty(obj,p)
            b = isempty(obj.AlignedWidth) && strcmp(p,'PaddingWidth');
        end

        % Display header
        function header = getHeader(obj)
            if ~isscalar(obj)
                header = getHeader@matlab.mixin.CustomDisplay(obj);
            else
                headerStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
                headerStr = ['  ', headerStr,' with field values:'];    % We will lead with the field names
                header = sprintf('%s\n',headerStr);
            end
        end
        
        % Returns property groups for display; fields, normal, read-only, and constant
        function groups = getPropertyGroups(obj)
       
            if ~isscalar(obj)
                groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
                return;
            end
        
            propName  = properties(obj);
            mc = meta.class.fromName(class(obj)); % Get metaclass
            allProperties = {mc.PropertyList.Name};
            % If defined explicitly, consider the order of properties as per 'CustomPropList'
            if any(strcmpi(allProperties,'CustomPropList'))
                propName = obj.CustomPropList;
            end
            propVal = cell(numel(propName),1);
            activeIdx = true(size(propName));
            nonFieldIdx = false(size(propName));  % Assume all are fields initially
            readOnlyIdx = false(size(propName));
            constIdx = false(size(propName));
            mc = meta.class.fromName(class(obj)); % Get metaclass
            % For each property...
            for n = 1:numel(propName)
                % Determine in property is read-only or constant from attributes
                id = cellfun(@(x)strcmp(x,propName{n}),{mc.PropertyList.Name});
                % If a read-only property then display it in a separate property group
                readOnlyIdx(n) = strcmp(mc.PropertyList(id).SetAccess,'private') && strcmp(mc.PropertyList(id).GetAccess,'public');
                % If a constant property then display it in a separate property group
                constIdx(n) = mc.PropertyList(id).Constant;  
                        
                % Don't use any active/inactive signalling here
                % if isInactiveProperty(obj,propName{n}) 
                %       % If it is inactive then do not add it to the property list
                %       activeIdx(n) = false;
                %  elseif isUndefinedProperty(obj,propName{n})
                %    % If it is undefined then set the value to an empty for disp()
                %        propVal{n} = [];
                %  else
                if isInactiveProperty(obj,propName{n}) 
                %       % If it is inactive then do not add it to the property list
                       activeIdx(n) = false;
                else
                % Capture the property value required and categorise the fields for display 
                if readOnlyIdx(n)
                    % Read-only property
                    propVal{n} = obj.(propName{n});         % Read-only, then get value
                else
                    % Writeable 
                    if isa(obj.(propName{n}),'BitField')
                        % Writeable BitField, treated as field by default
                        propVal{n} = obj.(propName{n}).Value;   % Assume that there is a 'Value' defined on this object
                    else
                        % Writable other, treat as field if another MessageFormat
                        nonFieldIdx(n) = ~isa(obj.(propName{n}),'MessageFormat');
                        propVal{n} = obj.(propName{n});
                    end
                end
                end
                % end
            end
            % Create separate property lists
            readOnlyPropList = cell2struct(propVal(activeIdx&readOnlyIdx),propName(activeIdx&readOnlyIdx));
            normalPropList = cell2struct(propVal(~nonFieldIdx&activeIdx&~readOnlyIdx&~constIdx),propName(~nonFieldIdx&activeIdx&~readOnlyIdx&~constIdx));
            nonFieldPropList = cell2struct(propVal(nonFieldIdx&activeIdx&~readOnlyIdx&~constIdx),propName(nonFieldIdx&activeIdx&~readOnlyIdx&~constIdx));
            constPropList = cell2struct(propVal(activeIdx&constIdx),propName(activeIdx&constIdx));
            % Combine the groups
            groups = [matlab.mixin.util.PropertyGroup(normalPropList) ...
                      matlab.mixin.util.PropertyGroup(nonFieldPropList,'Writeable properties:') ...
                      matlab.mixin.util.PropertyGroup(readOnlyPropList,'Read-only properties:')];                  
            if ~isempty(fields(constPropList))
                groups = [groups matlab.mixin.util.PropertyGroup(constPropList,'Constant properties:')];
            end
        end
            
    end
        
end

% Local functions

% Function is outside of width getter/class to get chain-through working through any 
% hierarchical message linking
function w = calcWidth(obj)
    [propName] = getActiveFieldProperties(obj,1);   % Do not include the local padding field at all
    w = 0;
    for n = 1:numel(propName)
        w = w+obj.(propName{n}).Width;
    end
end


