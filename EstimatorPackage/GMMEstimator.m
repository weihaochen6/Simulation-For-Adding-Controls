%% GMM estimator for time series
classdef GMMEstimator<handle
    %% Private Inputs Properties
    properties(Access = private)
        yData = []
        xData = []
        additional_controlData = []
        initialIVData = []
        % zData is a T x NumOfIV Matrix
        additionalIVData = []
    end

    %% Dependent Public Inputs Properties
    properties (Dependent, Access = public)
        % yData is T x 1 vector
        YData
        % xData is a T x K matrix
        XData
        Additional_controlData
        InitialIVData
        AdditionalIVData
    end

    %% Estimation Specification Properties
    properties(Access = public)
        % NoConstant is the option to control the constant or not
        NoConstant = true

        % The Significance Level when computing confidence intervals
        SignificanceLevel = 0.05

        % ErrorType specifies the error type
        ErrorType = 'iid'

        UseLongMoments = true
        InitialIV_AddinLongEquation = true
        AdditionalIV_AddinLongEquation = true
       
        UseShortMoments = false
        AdditionalIV_AddinShortEquation = false
        InitialIV_AddinShortEquation = false
    end

    %% OUTPUTS Properties
    properties(Access = public)
        FirstStepEstimates
        FirstStepVarCovMatrix
        FirstStepStandardErrors
        FirstSteptRatios
        FirstStepConfidenceInterval

        SecondStepEstimates
        SecondStepVarCovMatrix
        SecondStepStandardErrors
        SecondSteptRatios
        SecondStepConfidenceInterval
    end

    %% Basic Data Properties
    properties(GetAccess = public, SetAccess = private)
        % T is the number of time periods in the observed data
        T = 0
        % K is the number of x controls
        K = 0
        NumMomentsInShort = 0
        NumMomentsInLong = 0
        NumParInShort = 0
        NumParInLong = 0
        TotalNumOfMoments = 0
    end

    %% String Constants Properties
    properties (Constant,Access = private)
        % For Error Type
        iidString = 'iid'
        heteroskedasticityString = 'heteroskedasticity'
        AutocorrelatedString = 'autocorrelated'
    end

    %% HAC Estimator
    properties(GetAccess = public, SetAccess = private)
        HACEstimator = []
    end

    %% Some Intermediate Calculation Results
    properties(GetAccess = public, SetAccess = private)
        ShortControlMatrix= []
        LongControlMatrix= []
        ShortIVMatrix= []
        LongIVMatrix= []
        TotalIVMatrix= []

        ShortGradientOfMoments= []
        ShortMomentsAtZero= []
        LongGradientOfMoments= []
        LongMomentsAtZero = []
        TotalGradientOfMoment= []
        TotalMomentsAtZero= []
        FirstStepVarianceOfMoment= []
        SecondStepVarianceOfMoment= []

        FirstStepWeightingMatrix= []
        SecondStepWeightingMatrix= []

        FirstStepSigma2Hat= []
        SecondStepSigma2Hat= []
    end

    %% Constructor method
    methods(Access = public)
        function obj = GMMEstimator(inputStruct)
            % Allow create empty instance
            if nargin == 0
                return;
            end

            % Validate that the input argument is a structure
            if ~isstruct(inputStruct)
                error('GMMEstimator:Constructor:InputNotStruct', 'Input must be a structure.');
            end

            % Proceed if input is a structure
            propNames = fieldnames(inputStruct);
            for i = 1:length(propNames)
                propName = propNames{i};
                if isprop(obj, propName)
                    obj.(propName) = inputStruct.(propName);
                    %   else
                    % warning('GMMEstimator:Constructor:UnknownProperty', 'Property %s does not exist in the class.', propName);
                end
            end
        end
    end

    %%  set and get Methods for dependent properties
    methods
        %% For Inputs
        function set.YData(obj,YData)
            if ismatrix(YData) && size(YData,2) == 1
                obj.T = size(YData,1);
                obj.yData = YData;
            else
                error('YData must be a T x 1 vector.');
            end
        end

        function value = get.YData(obj)
            value = obj.yData;
        end

        function set.XData(obj,Data)
            if ismatrix(Data) && size(Data,1) == obj.T
                obj.K = size(Data,2);
                obj.xData = Data;
            else
                error('XData must be a T x K matrix and match the dimensions of YData.');
            end
        end

        function value = get.XData(obj)
            value = obj.xData;
        end

        function set.Additional_controlData(obj,Data)
            if ismatrix(Data) && size(Data,1) == obj.T
                obj.additional_controlData = Data;
            else
                error('XData must be a T x K matrix and match the dimensions of YData.');
            end
        end

        function value = get.Additional_controlData(obj)
            value = obj.additionalIVData;
        end

        function set.AdditionalIVData(obj,Data)
            if ismatrix(Data) && size(Data,1) == obj.T
                obj.additionalIVData = Data;
            else
                error('IV Data must be a T x NumOfIV matrix and match the time dimensions of YData.');
            end
        end

        function value = get.AdditionalIVData(obj)
            value = obj.additionalIVData;
        end

        function set.InitialIVData(obj,Data)
            if ismatrix(Data) && size(Data,1) == obj.T
                obj.initialIVData = Data;
            else
                error('InitialIVData must be a T x NumOfIV matrix and match the time dimensions of YData.');
            end
        end

        function value = get.InitialIVData(obj)
            value = obj.initialIVData;
        end
    end

    %% Main Public Methods
    methods(Access = public)
        %% First Step Estimation Results
        function FirstStepEstimates = getFirstStepEstimates(obj)
            if isempty(obj.yData) || isempty(obj.xData)
                error('Data is not imported');
            end

            % Set The Default Initial IV As x itself
            if isempty(obj.initialIVData)
                obj.initialIVData = obj.xData;
            end

            if obj.NoConstant
                ConstantVector = [];
            else
                ConstantVector = ones(obj.T,1);
            end


            % Long Equation
            obj.LongControlMatrix = [ConstantVector obj.xData obj.additional_controlData];
            obj.LongIVMatrix = [ConstantVector obj.additional_controlData];

            if obj.InitialIV_AddinLongEquation
                obj.LongIVMatrix = [obj.LongIVMatrix obj.initialIVData];
            end

            if obj.AdditionalIV_AddinLongEquation
                obj.LongIVMatrix = [obj.LongIVMatrix obj.additionalIVData];
            end

            if obj.UseLongMoments
                obj.NumParInLong = size(obj.LongControlMatrix,2);
                obj.NumMomentsInLong = size(obj.LongIVMatrix,2);
                obj.LongGradientOfMoments = (obj.LongIVMatrix'*obj.LongControlMatrix)/obj.T;
                obj.LongMomentsAtZero = (obj.LongIVMatrix'*obj.yData)/obj.T;
            end

            % Short Equation
            if obj.UseShortMoments
                obj.ShortControlMatrix = [ConstantVector obj.xData];
                obj.NumParInShort = size(obj.ShortControlMatrix,2);
                obj.NumMomentsInShort = size(obj.ShortIVMatrix,2);
                obj.ShortGradientOfMoments =obj.ShortIVMatrix'*obj.ShortControlMatrix/obj.T;
                obj.ShortMomentsAtZero = obj.ShortIVMatrix'*obj.yData/obj.T;

                if obj.InitialIV_AddinShortEquation
                    obj.ShortIVMatrix = [ConstantVector obj.initialIVData];
                end

                if obj.AdditionalIV_AddinShortEquation
                    obj.ShortIVMatrix = [obj.ShortIVMatrix obj.additionalIVData];
                end
            end

            if obj.UseLongMoments && ~obj.UseShortMoments
                obj.TotalMomentsAtZero = [obj.LongMomentsAtZero];
                obj.TotalGradientOfMoment = [obj.LongGradientOfMoments];
                obj.TotalIVMatrix = [obj.LongIVMatrix];
                obj.TotalNumOfMoments = obj.NumMomentsInLong;
            elseif ~obj.UseLongMoments && obj.UseShortMoments  
                obj.TotalMomentsAtZero = [obj.ShortMomentsAtZero];
                obj.TotalGradientOfMoment = [obj.ShortGradientOfMoments];
                obj.TotalIVMatrix = [obj.ShortIVMatrix];
                obj.TotalNumOfMoments = obj.NumMomentsInShort;
            elseif ~obj.UseLongMoments && ~obj.UseShortMoments  
                obj.TotalMomentsAtZero = [obj.ShortMomentsAtZero;obj.LongMomentsAtZero];
                obj.TotalGradientOfMoment = [obj.ShortGradientOfMoments zeros(obj.NumMomentsInShort, obj.NumParInLong-obj.NumParInShort); obj.LongGradientOfMoments];
                obj.TotalIVMatrix = [obj.ShortIVMatrix obj.LongIVMatrix];
                obj.TotalNumOfMoments = obj.NumMomentsInShort + obj.NumMomentsInLong;
            else
                error('You have to use one of moments from short equation or long equations.')
            end

            obj.FirstStepWeightingMatrix = eye(obj.TotalNumOfMoments)/(obj.TotalIVMatrix'*obj.TotalIVMatrix/obj.T);
            G = obj.TotalGradientOfMoment;
            denominator = (G')*obj.FirstStepWeightingMatrix*G;
            numerator = (G')*obj.FirstStepWeightingMatrix*(obj.TotalMomentsAtZero);
            obj.FirstStepEstimates = denominator\numerator;
            FirstStepEstimates = obj.FirstStepEstimates;
        end

        function FirstStepVarCovMatrix = getFirstStepVarCovMatrix(obj)
            if isempty(obj.FirstStepEstimates)
                getFirstStepEstimates(obj);
            end
            ShortResidual = [];
            if obj.UseShortMoments
                ShortEstimates = obj.FirstStepEstimates(1:obj.K);
                ShortResidual = obj.yData - obj.ShortControlMatrix*ShortEstimates;
            end

            LongResidual = [];
            if obj.UseLongMoments
                LongEstimates = obj.FirstStepEstimates;
                LongResidual = obj.yData - obj.LongControlMatrix*LongEstimates;
            end

            G = obj.TotalGradientOfMoment;
            W = obj.FirstStepWeightingMatrix;

            switch obj.ErrorType
                case obj.iidString
                    obj.FirstStepSigma2Hat = (LongResidual)'*(LongResidual)/(obj.T);
                    obj.FirstStepVarianceOfMoment =  obj.FirstStepSigma2Hat*eye(obj.TotalNumOfMoments)/(obj.TotalIVMatrix'*obj.TotalIVMatrix/obj.T);
                    obj.FirstStepVarCovMatrix = (1/obj.T)*obj.FirstStepSigma2Hat*eye(obj.NumParInLong)/(((G')*W*G));

                case obj.heteroskedasticityString
                    obj.FirstStepVarianceOfMoment = getVarianceOfMomentUnderHeteroskedasiticity(obj,ShortResidual,LongResidual);
                    Omega = obj.FirstStepVarianceOfMoment;
                    obj.FirstStepVarCovMatrix = (1/obj.T)*((G'*W*G)\((G')*W*Omega*(W')*G))/((G')*W*G);

                case obj.AutocorrelatedString
                    SampleMoments = obj.TotalIVMatrix.*ShortResidual;
                    obj.HACEstimator = NeweyWestHACEstimator(SampleMoments);
                    obj.FirstStepVarianceOfMoment = obj.HACEstimator.getLongRunVarianceEstimate();
                    Omega = obj.FirstStepVarianceOfMoment;
                    obj.FirstStepVarCovMatrix = (1/obj.T)*((G'*W*G)\((G')*W*Omega*(W')*G))/((G')*W*G);
            end

            FirstStepVarCovMatrix = obj.FirstStepVarCovMatrix;
        end

        function FirstStepStandardErrors = getFirstStepStandardErrors(obj)
            if isempty(obj.FirstStepVarCovMatrix)
                getFirstStepVarCovMatrix(obj);
            end

            obj.FirstStepStandardErrors = sqrt(diag(obj.FirstStepVarCovMatrix));
            FirstStepStandardErrors = obj.FirstStepStandardErrors;
        end

        function FirstSteptRatios = getFirstSteptRatios(obj)
            if isempty(obj.FirstStepStandardErrors)
                getFirstStepStandardErrors(obj);
            end

            obj.FirstSteptRatios = obj.FirstStepEstimates./obj.FirstStepStandardErrors;
            FirstSteptRatios = obj.FirstSteptRatios;
        end

        function ci = getFirstStepConfidenceInterval(obj)
            if isempty(obj.FirstStepStandardErrors)
                getFirstStepStandardErrors(obj);
            end
            c = norminv(1-obj.SignificanceLevel/2);
            beta = obj.FirstStepEstimates;
            se = obj.FirstStepStandardErrors;

            lb =  beta - c*se/sqrt(obj.T);
            ub =  beta + c*se/sqrt(obj.T);
            ci = [lb, ub];
            obj.FirstStepConfidenceInterval = ci;
        end

        %% First Step Auxiliary Functions
        function TestStatistic = getFirstStepTestStatistic(obj,k,Beta0)
            if isempty(obj.FirstStepStandardErrors)
                getFirstStepStandardErrors(obj);
            end
            if nargin < 2
                k = 1;
            end
            if nargin < 3
                Beta0 = 0;
            end
            beta = obj.FirstStepEstimates;
            se = obj.FirstStepStandardErrors(k);
            TestStatistic = (beta(k) - Beta0)/se;
        end

        function TestConclusion = FirstStepTestOfCoefficient(obj,k,Beta0)
            % TestConclusion = 1 means rejection
            if isempty(obj.FirstStepStandardErrors)
                getFirstStepStandardErrors(obj);
            end
            if nargin < 2
                k = 1;
            end
            if nargin < 3
                Beta0 = 0;
            end

            StepTestStatistic = getFirstStepTestStatistic(obj,k,Beta0);

            CriticalValue = norminv(1-obj.SignificanceLevel/2);
            if(abs(StepTestStatistic)>CriticalValue)
                TestConclusion = 1;
            else
                TestConclusion = 0;
            end
        end

        function ESC = getFirstStepESC(obj)
            if isempty(obj.FirstStepVarCovMatrix)
                getFirstStepVarCovMatrix(obj);
            end
            V = obj.FirstStepVarCovMatrix(1,1)*obj.T;
            c = obj.TotalNumOfMoments;
            ESC = log(V) + (c-obj.K)*log(sqrt(obj.T))/sqrt(obj.T);
        end

        %% Second Step Estimation Results
        function SecondStepEstimates = getSecondStepEstimates(obj)
            if isempty(obj.FirstStepEstimates)
                obj.getFirstStepEstimates();
                obj.getFirstStepVarCovMatrix();
            end
            obj.SecondStepWeightingMatrix = eye(obj.TotalNumOfMoments)/obj.FirstStepVarianceOfMoment;
            
            
            G = obj.TotalGradientOfMoment;
            W = obj.SecondStepWeightingMatrix;

            denominator = (G')*W*G;
            numerator = (G')*W*obj.TotalMomentsAtZero;
            obj.SecondStepEstimates = denominator\numerator;
            SecondStepEstimates = obj.SecondStepEstimates;
            
        end

        function SecondStepVarCovMatrix = getSecondStepVarCovMatrix(obj)
            if isempty(obj.SecondStepEstimates)
                getSecondStepEstimates(obj);
            end

            ShortResidual = [];
            if obj.UseShortMoments
                ShortEstimates = obj.SecondStepEstimates(1:obj.K);
                ShortResidual = obj.yData - obj.ShortControlMatrix*ShortEstimates;
            end

            LongResidual = [];
            if obj.UseLongMoments
                LongEstimates = obj.SecondStepEstimates;
                LongResidual = obj.yData - obj.LongControlMatrix*LongEstimates;
            end

            G = obj.TotalGradientOfMoment;
            W = obj.SecondStepWeightingMatrix;

            switch obj.ErrorType
                case obj.iidString
                    obj.SecondStepSigma2Hat = (LongResidual)'*(LongResidual)/(obj.T);
                    obj.SecondStepVarianceOfMoment = obj.SecondStepSigma2Hat*eye(obj.TotalNumOfMoments)/(obj.TotalIVMatrix'*obj.TotalIVMatrix/obj.T);
                    obj.SecondStepVarCovMatrix = (1/obj.T)*obj.SecondStepSigma2Hat*eye(obj.NumParInLong)/(((G')*W*G));

                case obj.heteroskedasticityString
                    obj.SecondStepVarianceOfMoment = getVarianceOfMomentUnderHeteroskedasiticity(obj,ShortResidual,LongResidual);
                    Omega = obj.SecondStepVarianceOfMoment;
                    obj.SecondStepVarCovMatrix = (1/obj.T)*((G'*W*G)\((G')*W*Omega*(W')*G))/((G')*W*G);

                case obj.AutocorrelatedString
                    SampleMoments = obj.TotalIVMatrix.*Residual;
                    obj.HACEstimator = NeweyWestHACEstimator(SampleMoments);
                    obj.SecondStepVarianceOfMoment = obj.HACEstimator.getLongRunVarianceEstimate();
                    Omega = obj.SecondStepVarianceOfMoment;
                    obj.SecondStepVarCovMatrix = (1/obj.T)*((G'*W*G)\((G')*W*Omega*(W')*G))/((G')*W*G);
            end

            SecondStepVarCovMatrix = obj.SecondStepVarCovMatrix;
        end

        function SecondStepStandardErrors = getSecondStepStandardErrors(obj)
            if isempty(obj.SecondStepVarCovMatrix)
                getSecondStepVarCovMatrix(obj);
            end

            obj.SecondStepStandardErrors = sqrt(diag(obj.SecondStepVarCovMatrix));
            SecondStepStandardErrors = obj.SecondStepStandardErrors;
        end

        function SecondSteptRatios = getSecondSteptRatios(obj)
            if isempty(obj.SecondStepStandardErrors)
                getSecondStepStandardErrors(obj);
            end
            obj.SecondSteptRatios = obj.SecondStepEstimates./obj.SecondStepStandardErrors;
            SecondSteptRatios = obj.SecondSteptRatios;
        end

        function ci = getSecondStepConfidenceInterval(obj)
            if isempty(obj.SecondStepStandardErrors)
                getSecondStepStandardErrors(obj);
            end
            c = norminv(1-obj.SignificanceLevel/2);
            beta = obj.SecondStepEstimates;
            se = obj.SecondStepStandardErrors;

            lb =  beta - c*se/sqrt(obj.T);
            ub =  beta + c*se/sqrt(obj.T);
            ci =  [lb, ub];
            obj.SecondStepConfidenceInterval = ci;
        end

        %% Second Step Auxiliary Functions
        function TestStatistic = getSecondStepTestStatistic(obj,k,Beta0)
            if isempty(obj.SecondStepStandardErrors)
                getSecondStepStandardErrors(obj);
            end
            if nargin < 2
                k = 1;
            end
            if nargin < 3
                Beta0 = 0;
            end
            beta = obj.SecondStepEstimates;
            se = obj.SecondStepStandardErrors(k);
            TestStatistic = (beta(k) - Beta0)/se;
        end

        function TestConclusion = SecondStepTestOfCoefficient(obj,k,Beta0)
            % TestConclusion = 1 means rejection
            if isempty(obj.SecondStepStandardErrors)
                getSecondStepStandardErrors(obj);
            end
            if nargin < 2
                k = 1;
            end
            if nargin < 3
                Beta0 = 0;
            end
            TestStatistic = getSecondStepTestStatistic(obj,k,Beta0);

            CriticalValue = norminv(1-obj.SignificanceLevel/2);
            if(abs(TestStatistic)>CriticalValue)
                TestConclusion = 1;
            else
                TestConclusion = 0;
            end
        end

        function ESC = getSecondStepESC(obj)
            if isempty(obj.SecondStepVarCovMatrix)
                getSecondStepVarCovMatrix(obj);
            end
            V = obj.SecondStepVarCovMatrix(1,1)*obj.T;
            c = obj.TotalNumOfMoments;
            ESC = log(V) + (c-obj.K)*log(sqrt(obj.T))/sqrt(obj.T);
        end

    end

    %% Private Methods For Calculation
    methods(Access = private)
        function VarianceOfMoment = getVarianceOfMomentUnderHeteroskedasiticity(obj,ShortResiduals,LongResiduals)
            LongIndividualMoments = LongResiduals.*obj.LongIVMatrix;
            ShortIndividualMoments = [];
            if obj.UseShortMoments
                ShortIndividualMoments  = ShortResiduals.*obj.ShortIVMatrix;
            end
            TotalIndividualMoments = [ShortIndividualMoments LongIndividualMoments];
            VarianceOfMoment = (TotalIndividualMoments'*TotalIndividualMoments)/obj.T;
        end
    end

    %% Methods For Simulations
    methods(Access = public)
        function obj = FinalizeEstimation(obj)
            % Clear The Data
            obj.yData = [];
            obj.xData = [];
            obj.additionalIVData = [];
            obj.additional_controlData = [];
            obj.initialIVData = [];
            % Clear Some Intermediate Calculation Results
            obj.ShortControlMatrix= [];
            obj.LongControlMatrix= [];
            obj.ShortIVMatrix= [];
            obj.LongIVMatrix= [];
            obj.TotalIVMatrix= [];
            obj.TotalIVMatrix= [];
        end
    end
end
