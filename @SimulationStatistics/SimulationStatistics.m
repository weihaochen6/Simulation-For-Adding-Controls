classdef SimulationStatistics < handle
    %% Inputs Properties
    properties(Access = public)
        EstimationResults = []
        TrueValues

        % Round the results
        RoundOption = 4
    end
    
    %% OUTPUTS Properties
    properties(GetAccess = public, SetAccess = public)
        Mean
        MeanBias
        Median
        MedianBias
        SampleVariance
        MSE
        RMSE
        % Median Absolute Error 
        MAE
    end

    %% Constructor method
    methods(Access = public)
        function obj = SimulationStatistics(inputStruct)
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
    
    %% Main Public Methods
    methods(Access = public)
        function meanValue = getMean(obj, fieldname)
            % Use extractValues to unify value extraction
            values = extractValues(obj, fieldname);
            obj.Mean = mean(values, 2); % Compute the mean along columns
            meanValue = obj.roundResults(obj.Mean);
        end

        function meanBiasValue = getMeanBias(obj, fieldname)
            if ~isempty(obj.TrueValues)
                meanValue = getMean(obj, fieldname);
                obj.MeanBias = meanValue - obj.TrueValues; % Assign mean value to Mean property
            else
                obj.MeanBias = [];
            end
            meanBiasValue =  obj.roundResults(obj.MeanBias);
        end

        function MedianValue = getMedian(obj, fieldname)
            % Extract values using the extractValues function
            values = extractValues(obj, fieldname);
            obj.Median  = median(values, 2); % Compute the median along columns
            MedianValue = obj.roundResults(obj.Median); % Assign median value to Median property
        end

        function MedianBiasValue = getMedianBias(obj,fieldname)
            MedianValue = getMedian(obj, fieldname);
            if ~isempty(obj.TrueValues)
                obj.MedianBias = MedianValue - obj.TrueValues; % Assign median bias value to MedianBias property
            else
                obj.MedianBias = []; % Set MedianBias to empty if TrueValues is not set
            end
            MedianBiasValue = obj.roundResults(obj.MedianBias);
        end

        function SampleVariance = getSampleVariance(obj, fieldname)
            % Extract values using the extractValues function
            values = extractValues(obj, fieldname);
            obj.SampleVariance = var(values, 0, 2); % Compute the sample variance along columns
             SampleVariance = obj.roundResults(obj.SampleVariance); % Assign sample variance value to SampleVariance property
        end

         function MAEValue = getMAE(obj, fieldname)
            if ~isempty(obj.TrueValues)
                estimates = obj.extractValues(fieldname);
                deviations = abs(bsxfun(@minus, estimates, obj.TrueValues));
                obj.MAE = median(deviations, 2);
            else
                obj.MAE = [];
            end
            MAEValue = obj.roundResults(obj.MAE);
        end


        function MSEValue = getMSE(obj, fieldname)
            if ~isempty(obj.TrueValues)
                estimates = obj.extractValues(fieldname);
                deviations = (estimates - obj.TrueValues).^2;
                obj.MSE = mean(deviations, 2);
            else
                obj.MSE = [];
            end
            MSEValue = obj.roundResults(obj.MSE);
        end

         function RMSEValue = getRMSE(obj, fieldname)
            if isempty(obj.MSE)
                getMSE(obj, fieldname);
                obj.RMSE = sqrt(obj.MSE);
            else
                 obj.RMSE = [];
            end
            RMSEValue = obj.roundResults(obj.RMSE);
        end
    end

    %% Private Methods
    methods(Access = private)
        function values = extractValues(obj, fieldname)
            % Utility function to extract and concatenate field values
            if isempty(obj.EstimationResults) || ~isprop(obj.EstimationResults(1), fieldname)
                error([class(obj) ':extractValues:InvalidFieldName'], ...
                    ['Estimation Results Not Imported or Field ' fieldname ' does not exist in the objects of EstimationResults.']);
            end
            values = arrayfun(@(x) x.(fieldname), obj.EstimationResults, 'UniformOutput', false);
            values = [values{:}];
            if isempty(values)
                error([class(obj) ':extractValues:EmptyFieldValues'], ...
                    'All objects have empty values for the specified field.');
            end
        end

        function result = roundResults(obj, value)
            if isempty(obj.RoundOption)
                result = value; % Do not round
            elseif isinteger(obj.RoundOption)
                result = round(value, obj.RoundOption); % Round to specified digits
            else
                error([class(obj) ':roundResults:InvalidNumRound'], ...
                    'Num_Round must be an integer or empty.');
            end
        end
    end

    %% Methods For Simulation
    methods(Access = public)
        function FinalizeCalculation(obj)
            % Set EstimationResults property to empty
            obj.EstimationResults = [];
        end
    end
end
