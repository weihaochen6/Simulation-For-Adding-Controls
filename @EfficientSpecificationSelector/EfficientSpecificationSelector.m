classdef EfficientSpecificationSelector<handle
    %% Private Inputs Properties
    properties(Access = private)
        % ESCArray is an array storing ESC
        escArray = []
    end

    %% Dependent Public Inputs Properties
    properties(Dependent, Access = public)
        % ESCArray is an array storing ESC
        ESCArray
    end

    %% OUTPUTS properties
    properties (GetAccess = public,SetAccess = private)
        ESC_Minimum
        LinearIndex
    end

    %% Constructor Methods
    methods(Access = public)
        function obj = EfficientSpecificationSelector(inputStruct)
            % Allow create empty instance
            if nargin == 0
                return;
            end

            % Validate that the input argument is a structure
            if ~isstruct(inputStruct)
                error([class(obj) ':Constructor:InputNotStruct'], 'Input must be a structure.');
            end

            % Proceed if input is a structure
            propNames = fieldnames(inputStruct);
            for i = 1:length(propNames)
                propName = propNames{i};
                if isprop(obj, propName)
                    obj.(propName) = inputStruct.(propName);
                end
            end
        end
    end
    %%  set and get Methods for dependent properties
    methods
        function set.ESCArray(obj,ESCArray)
            obj.escArray = ESCArray;
            obj.ESC_Minimum = [];
            obj.LinearIndex = [];
        end

        function value = get.ESCArray(obj)
            value = obj.escArray;
        end
    end

    %% Main Public Methods
    methods(Access = public)
        function ESC_Minimum = getESC_Minimum(obj)
            [ESC_Minimum, obj.LinearIndex] = min(obj.ESCArray, [], 'all', 'linear');
        end

        function varargout = getIndexOfMinimum(obj)
            % Ensure the index of the minimum value is calculated
            if isempty(obj.LinearIndex)
                obj.getESC_Minimum();
            end

            % Get the number of dimensions of ESCArray
            sz = size(obj.ESCArray);
            nDims = length(sz); % Get the number of dimensions

            % Check if the number of outputs matches the number of dimensions
            if nargout ~= nDims
                error('The number of output arguments must match the number of dimensions in ESCArray.');
            end

            % Use a cell array to receive the multiple outputs of ind2sub
            indicesCell = cell(1, nDims);
            [indicesCell{:}] = ind2sub(sz, obj.LinearIndex);

            % Prepare varargout output, determine the number of outputs based on nargout
            varargout = cell(1, nargout); % At this point, nargout should be equal to nDims

            % Assign the outputs of ind2sub to varargout
            for k = 1:nargout
                varargout{k} = indicesCell{k};
            end
        end

    end
end