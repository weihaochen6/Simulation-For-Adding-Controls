%% GMM estimator for time series
classdef TSLSEstimator<handle
    %% Inputs Properties
    properties (Access = public)
        yData = []
        yLagNum = 0
        xData = []

        % Data For IV
        zData = []
        NoConstant = true
        ErrorType = 'iid'
    end

    %% OUTPUTS Properties
    properties(Access = public)
        FirstStageEstimates

        Estimates
        VarCovMatrix
        StandardErrors
        tRatios
    end

    %% Basic Data Properties
    properties(GetAccess = public, SetAccess = private)
        T    % This will store the number of time periods
        K    % This will store the number of x regressors
        NumOfIV
        NumOfRegressors
    end

    %% Properties for Clarifications of The Code
    properties(Access = public)
        % properties(Access = private)
        % Two Values: 'No X regressors' or 'With X regressors'
        WithOrNoXregressorsCases = 'With X regressors'
    end

    %% String Constants Properties
    properties (Constant,Access = private)
        WithXregressorsStringConstant = 'With X regressors'
        NoXregressorsStringConstant = 'No X regressors'
    end

    properties(GetAccess = public, SetAccess = private)
        ValidNumOfT % The real number of time periods in the model

        % GetSubVectors
        yVector
        yLagVector
        XVector
        ZVector

        WMatrix
        PredictedWMatrix
        IVMatrix

        Sigma2Hat
    end

    %% Constructor method
    methods(Access = public)
        function obj = TSLSEstimator(YData, XData, ZData)
            % Allow create empty instance
            if nargin == 0
                return;
            end
            % Handle YData
            if(nargin>=1 && ~isempty(YData))
                if obj.SetyData(YData)
                    obj.yData = YData;
                end
            end
            % Handle XData
            if nargin >= 2 && ~isempty(XData)
                if obj.SetxData(XData)
                    obj.xData = XData;
                    obj.WithOrNoXregressorsCases = obj.WithXregressorsStringConstant;
                end
            end

            % Handle ZData
            if nargin >=3 && ~isempty(ZData)
                if obj.SetzData(ZData)
                    obj.zData = ZData;
                end
            end
        end
    end

    %% Public Set Methods
    methods
        function set.yData(obj,YData)
            if obj.SetyData(YData)
                obj.yData = YData;
            end
        end

        function set.xData(obj,XData)
            if obj.SetxData(XData)
                obj.xData = XData;
            end
        end
        function set.zData(obj,ZData)
            if obj.SetzData(ZData)
                obj.zData = ZData;
            end
        end
    end

    %% Private Set Methods
    methods(Access = private)
        function isValid = SetyData(obj, YData)
            if ismatrix(YData) && size(YData,2) == 1
                isValid = true;
                obj.T = size(YData,1);
            else
                error('YData must be a T x 1 vector.');
            end
        end

        function isValid = SetxData(obj, XData)
            if ismatrix(XData) && size(XData,1) == obj.T
                isValid = true;
                obj.K = size(XData,2);
            else
                error('XData must be a T x K matrix and match the dimensions of YData.');
            end
        end

        function isValid = SetzData(obj,ZData)
            if ismatrix(ZData) && size(ZData,1) == obj.T
                isValid = true;
                obj.NumOfIV = size(ZData,2);
            else
                error('IV Data must be a T x NumOfIV matrix and match the time dimensions of YData.');
            end
        end
    end

    %% Main Public Methods
    methods(Access = public)
        %% First Step Results
        function Estimates = getEstimates(obj)
            obj.ValidNumOfT = obj.T - obj.yLagNum;
            GetSubVectors(obj);
            ConstructWMatrix(obj);
            ConstructIVMatrix(obj);
            Y = obj.yVector;
            
            obj.FirstStageEstimates = obj.WMatrix\obj.IVMatrix;
            obj.PredictedWMatrix = obj.WMatrix*obj.FirstStageEstimates;
            obj.Estimates = (obj.PredictedWMatrix'*obj.WMatrix)\(obj.PredictedWMatrix'*Y);
            
            obj.NumOfRegressors = length(obj.Estimates);
            Estimates = obj.Estimates;
        end

        function VarCovMatrix = getVarCovMatrix(obj)
            if isempty(obj.Estimates)
                getEstimates(obj);
            end
            Residual = obj.yVector - obj.WMatrix*obj.Estimates;
            obj.Sigma2Hat = (Residual)'*(Residual)/(obj.ValidNumOfT);

            obj.VarCovMatrix = obj.Sigma2Hat*eye(obj.NumOfRegressors)/(obj.PredictedWMatrix'*obj.PredictedWMatrix);
            VarCovMatrix = obj.VarCovMatrix;
        end
        function StandardErrors = getStandardErrors(obj)
            if isempty(obj.VarCovMatrix)
                getVarCovMatrix(obj);
            end

            obj.StandardErrors = sqrt(diag(obj.VarCovMatrix));
            StandardErrors = obj.StandardErrors;
        end

        function tRatios = gettRatios(obj)
            if isempty(obj.StandardErrors)
                getStandardErrors(obj);
            end
            obj.tRatios = obj.Estimates./obj.StandardErrors;
            tRatios = obj.FirstSteptRatios;
        end
    end
    %% Private Methods For Calculation
    methods(Access = private)
        function GetSubVectors(obj)
            obj.yVector = obj.yData(obj.yLagNum+1:end);
            if obj.yLagNum>=1
                obj.yLagVector = zeros(obj.ValidNumOfT,obj.yLagNum);
                for k = 1:obj.yLagNum
                    obj.yLagVector(:,k) = obj.yData(obj.yLagNum+1-k:end-k);
                end
            end
            obj.XVector = obj.xData(obj.yLagNum+1:end,:);
            obj.ZVector = obj.zData(obj.yLagNum+1:end,:);
        end

        function obj = ConstructWMatrix(obj)
            if(obj.NoConstant)
                obj.WMatrix = [obj.yLagVector,obj.XVector];
            else
                obj.WMatrix = [ones(obj.ValidNumOfT,1) obj.yLagVector obj.XVector];
            end
        end

        function obj = ConstructIVMatrix(obj)
            if(obj.NoConstant)
                obj.IVMatrix = [obj.yLagVector,obj.XVector obj.ZVector];
            else
                obj.IVMatrix = [ones(obj.ValidNumOfT,1) obj.yLagVector obj.XVector obj.ZVector];
            end
        end
    end
end