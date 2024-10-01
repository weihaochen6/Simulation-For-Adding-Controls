% This class is used to estimate the
classdef AdaptiveLassoSelector<handle
    %% Private Inputs Properties
    properties(Access = private)
        residuals = []
        zData = []
        candidateData = []
    end

    %% Dependent Public Inputs Properties
    properties (Dependent, Access = public)
        Residuals
        ZData
        CandidateData
    end

    %% Selection Setting Properties
    properties(Access = public)
        MSCType= 'BIC'
        LambdaStep = 0.01
        MaxLambda = 1
        Tau = 1
        True_z_Uncor
        True_u_Uncor
    end

    %% String Constants For Specifications
    properties(Constant,Access = private)
        BICStr = 'BIC'
        AICStr = 'AIC'
    end

    %% Basic Data Properties
    properties(GetAccess = public, SetAccess = private)
        T
        Num_z
        Num_Candidate
    end

    %% OUTPUTS properties
    properties (GetAccess = public,SetAccess = private)
        zUncorMinimizer %  K x Num_Candidate x NumOfLambdaValues
        ResidualUncorMinimizer % Num_Candidate x NumOfLambdaValues

        SelectedzUncor  % Num_Candidate x 1
        SelectedResidualUncor % Num_Candidate x 1

        % Identity Vector For Four Types 
        SelectedzUncor_ResidualUncor
        SelectedzUncor_ResidualCor
        SelectedzCor_ResidualUncor
        SelectedzCor_ResidualCor

        % Identity Vector For Those Should Be Added As IV and Control
        Selected_Candidate_IV
        Num_Selected_Candidate_IV
        % Identity Vector For Those Should Be Added As Control
        Selected_Candidate_Control
        Num_Selected_Candidate_Control

        SelectedLambdazUncor
        SelectedLambdaResidualUncor

        GMM_MSC_zUncor %  Num_Candidate x 1 vector
        GMM_MSC_ResidualUncor %  Num_Candidate x 1 vector
    end

    %% Private Calculation Properties
    properties (Access = private)
        LambdaRange
        Num_LambdaValues
        Expectation_Candi_Residual
        Expectation_Candi_z
        L2Norm_Expectation_Candi_z
        L2Norm_Expectation_Candi_Residual
        Penalty
    end

    %% Constructor Methods
    methods(Access = public)
        function obj = AdaptiveLassoSelector(inputStruct)
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
    %% Set and Get Methods
    % Public Set Methods
    methods
        function set.Residuals(obj,Residuals)
            if ~isvector(Residuals) || size(Residuals,1)<=1
                error('Residuals must be a T x 1 vector.');
            end

            obj.T = length(Residuals);
            obj.residuals = Residuals;
        end

        function value = get.Residuals(obj)
            value = obj.residuals;
        end

        function set.ZData(obj,zData)
            if ~isvector(zData) && ~ismatrix(zData)
                error('ControlData must be a vector or matrix.');
            end
            obj.Num_z = size(zData, 2);
            obj.zData = zData;
        end

        function value = get.ZData(obj)
            value = obj.zData;
        end

        function set.CandidateData(obj,data)
            if ~isvector(data) && ~ismatrix(data)
                error('IV must be a vector or matrix.');
            end
            obj.Num_Candidate = size(data,2);
            obj.candidateData = data;
        end

        function value = get.CandidateData(obj)
            value = obj.candidateData;
        end
    end

    %% Main Public Methods
    methods(Access = public)
        function obj = getSelectionResults(obj)
            obj.Expectation_Candi_Residual = obj.candidateData'*obj.residuals/obj.T;
            obj.Expectation_Candi_z = obj.candidateData'*obj.zData/obj.T;
            obj.L2Norm_Expectation_Candi_z = abs(obj.Expectation_Candi_z);
            obj.L2Norm_Expectation_Candi_Residual = abs(obj.Expectation_Candi_Residual);

            obj.LambdaRange = 0:obj.LambdaStep:obj.MaxLambda;
            obj.Num_LambdaValues = length(obj.LambdaRange);

            switch obj.MSCType
                case obj.BICStr
                    obj.Penalty = -log(obj.T);
                case obj.AICStr
                    obj.Penalty = -2;
                otherwise
                    error('wrong moment selection criterion input.');
            end

            obj.zUncorMinimizer = zeros(obj.Num_Candidate, obj.Num_z,obj.Num_LambdaValues);

            obj.ResidualUncorMinimizer = zeros(obj.Num_Candidate,obj.Num_LambdaValues);

            obj.GMM_MSC_zUncor = zeros(obj.Num_LambdaValues,1);
            obj.GMM_MSC_ResidualUncor = zeros(obj.Num_LambdaValues,1);

            for i = 1:obj.Num_LambdaValues
                lambda = obj.LambdaRange(i);
                [obj.zUncorMinimizer(:,:,i), obj.GMM_MSC_zUncor(i)] = Select_IVUncorrelatedOnes(obj,lambda);
                [obj.ResidualUncorMinimizer(:,i), obj.GMM_MSC_ResidualUncor(i)] = Select_ErrorUncorrelatedOnes(obj,lambda);
            end

            % getSelectedIdentities
            % For IV Uncorrelated Part
            [~,index] = min(obj.GMM_MSC_zUncor);
            obj.SelectedLambdazUncor = obj.LambdaRange(index);

            %obj.SelectedzUncor  = transpose(double(all(obj.zUncorMinimizer(:,:,index)==0)));
            obj.SelectedzUncor  = double(obj.zUncorMinimizer(:,:,index)==0);

            % For Error Uncorrelated Part
            [~,index] = min(obj.GMM_MSC_ResidualUncor);
            obj.SelectedLambdaResidualUncor = obj.LambdaRange(index);
            obj.SelectedResidualUncor  = double(obj.ResidualUncorMinimizer(:,index) == 0);

            %% get The Identity Vector For Each Types
            obj.SelectedzUncor_ResidualUncor =  obj.SelectedzUncor.*obj.SelectedResidualUncor;
            obj.SelectedzUncor_ResidualCor = obj.SelectedzUncor.*(1-obj.SelectedResidualUncor);
            obj.SelectedzCor_ResidualUncor = (1-obj.SelectedzUncor).*obj.SelectedResidualUncor;
            obj.SelectedzCor_ResidualCor = (1-obj.SelectedzUncor).*(1-obj.SelectedResidualUncor);

            obj.Selected_Candidate_IV = obj.SelectedzUncor_ResidualUncor + obj.SelectedzCor_ResidualUncor;
            obj.Num_Selected_Candidate_IV = sum(obj.Selected_Candidate_IV);
            obj.Selected_Candidate_Control = obj.SelectedzUncor_ResidualCor;
            obj.Num_Selected_Candidate_Control = sum(obj.Selected_Candidate_Control);
        end

        %% Auxiliary Public Functions 
        function value = isAllTrue_z_UncorSelected(obj)
            if all(obj.SelectedzUncor == obj.True_z_Uncor)
                value = 1;
            else
                value = 0;
            end
        end

        function value = isAllTrue_u_UncorSelected(obj)
            if all(obj.SelectedResidualUncor == obj.True_u_Uncor)
                value = 1;
            else
                value = 0;
            end
        end

        function value = isAllTrue_z_Uncor_True_u_UncorSelected(obj)
            if all(obj.SelectedzUncor == obj.True_z_Uncor) && all(obj.SelectedResidualUncor == obj.True_u_Uncor)
                value = 1;
            else 
                value = 0;
            end
        end

        function [RowNum, ColNum] = getDimensionsForEstimationResults(obj)
            RowNum = 2^obj.Num_Selected_Candidate_IV;
            ColNum = 2^obj.Num_Selected_Candidate_Control;
        end

        function SelectionCountMatrix = getSelectionCountMatrix(obj)
            SelectionCountMatrix = zeros(2^obj.Num_Selected_Candidate_IV,2^obj.Num_Selected_Candidate_Control);
            [rowIndex, colIndex] = getIndicesOfSelected_On_2_to_NumCandidate_SquareMatrix(obj);
            SelectionCountMatrix(rowIndex, colIndex) = 1;
        end

        function  SelectedControlData = getSelectedControlData(obj)
            SelectedControlData = obj.candidateData(:,logical(obj.Selected_Candidate_Control));
        end

        function SelectedIVData = getSelectedIVData(obj)
            SelectedIVData = obj.candidateData(:,logical(obj.Selected_Candidate_IV));
        end
    end

    %% private calculation methods
    methods(Access = private)
        function [ivUncorMinimizer, gmmMSCForIVUncor] = Select_IVUncorrelatedOnes(obj,lambda)
            lambda = lambda./((obj.L2Norm_Expectation_Candi_z).^(obj.Tau));
            ivUncorMinimizer = sign(obj.Expectation_Candi_z).*(max(obj.L2Norm_Expectation_Candi_z - lambda,0));
            df = sum(ivUncorMinimizer==0);
            gmmMSCForIVUncor = obj.T*((obj.Expectation_Candi_z - ivUncorMinimizer)'*(obj.Expectation_Candi_z - ivUncorMinimizer))+...
                obj.Penalty*df; 
        end

        function [rowIndex, colIndex] = getIndicesOfSelected_On_2_to_NumCandidate_SquareMatrix(obj)
            vec1 = obj.Selected_Candidate_IV;
            vec2 = obj.Selected_Candidate_Control;
            weights = 2.^(obj.Num_Candidate-1:-1:0);
            rowIndex = sum(weights*vec1)+1;
            colIndex = sum(weights*vec2)+1;
        end

        function [errorUncorMinimizer, gmmMSCForErrorUncor] = Select_ErrorUncorrelatedOnes(obj,lambda)
            lambda = lambda./((obj.L2Norm_Expectation_Candi_Residual).^(obj.Tau));
            errorUncorMinimizer = sign(obj.Expectation_Candi_Residual).*(max(obj.L2Norm_Expectation_Candi_Residual - lambda,0));% num_control x 1
            df = sum(errorUncorMinimizer==0);
            gmmMSCForErrorUncor = obj.T*((obj.Expectation_Candi_Residual - errorUncorMinimizer)'*(obj.Expectation_Candi_Residual - errorUncorMinimizer))+obj.Penalty*df;
        end
    end



    %% Methods For Simulations
    methods(Access = public)
        function FinalizeSelection(obj)
            obj.residuals = [];
            obj.zData = [];
            obj.candidateData = [];
        end
    end
end
