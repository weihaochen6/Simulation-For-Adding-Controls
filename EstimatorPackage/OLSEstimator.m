classdef OLSEstimator<handle
    %% Inputs Properties
    properties (Access = public)
        yData = []
        yLagNum = 0
        xData = []
        NoConstant = true
        ErrorType = 'iid' 
    end

    %% OUTPUTS Properties
    properties(Access = public)
        Estimates
        VarCovMatrix
        StandardErrors
        tRatios
    end

    %% Basic Data Properties
    properties(GetAccess = public, SetAccess = private)
        T    % This will store the number of time periods
        K    % This will store the number of x regressors
        NumOfRegressors
    end

    %% Properties for Clarifications of The Code
    properties(Access = public)
        % properties(Access = private)
        % Two Values: 'No X regressors' or 'With X regressors'
        WithOrNoXregressorsCases = 'No X regressors'
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
        WMatrix
        Sigma2Hat
    end

    %% Constructor method
    methods(Access = public)
        function obj = OLSEstimator(YData, XData)
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
    end


    %% Main Public Methods
    methods(Access = public)
        function Estimates = getEstimates(obj)
            obj.ValidNumOfT = obj.T - obj.yLagNum;
            GetSubVectors(obj);
            ConstructWMatrix(obj);
            Y = obj.yVector;
            obj.Estimates = obj.WMatrix\Y;
            obj.NumOfRegressors = length(obj.Estimates);
            Estimates = obj.Estimates;
        end

        function VarCovMatrix = getVarCovMatrix(obj)
            if isempty(obj.Estimates)
                getEstimates(obj);
            end
            Residual = obj.yVector - obj.WMatrix*obj.Estimates;
            obj.Sigma2Hat = (Residual)'*(Residual)/(obj.ValidNumOfT - obj.NumOfRegressors);
            obj.VarCovMatrix = obj.Sigma2Hat * eye*(obj.NumOfRegressors)/(obj.WMatrix'*obj.WMatrix);
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
            tRatios = obj.tRatios;
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
        end
        function obj = ConstructWMatrix(obj)
            if(obj.NoConstant)
                obj.WMatrix = [obj.yLagVector,obj.XVector];
            else
                obj.WMatrix = [ones(obj.ValidNumOfT,1) obj.yLagVector obj.XVector];
            end
        end
    end
end