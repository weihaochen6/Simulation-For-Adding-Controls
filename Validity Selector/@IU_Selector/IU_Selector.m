% This class is used to estimate the
classdef IU_Selector<handle
    %% Private Inputs Properties
    properties(Access = private)
        zData = []
        candidateData = []
    end

    %% Dependent Public Inputs Properties
    properties (Dependent, Access = public)
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
        SelectedzUncor  % Num_Candidate x 1

        % Identity Vector For Those Should Be Added As Control
        Num_Selected_Candidate_Control

        SelectedLambdazUncor
        GMM_MSC_zUncor %  Num_Candidate x 1 vector
    end

    %% Private Calculation Properties
    properties (Access = private)
        LambdaRange
        Num_LambdaValues
        Expectation_Candi_z
        L2Norm_Expectation_Candi_z
        Penalty
    end

    %% Constructor Methods
    methods(Access = public)
        function obj = IU_Selector(inputStruct)
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
            obj.Expectation_Candi_z = obj.candidateData'*obj.zData/obj.T;
            obj.L2Norm_Expectation_Candi_z = abs(obj.Expectation_Candi_z);

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


            obj.GMM_MSC_zUncor = zeros(obj.Num_LambdaValues,1);

            for i = 1:obj.Num_LambdaValues
                lambda = obj.LambdaRange(i);
                [obj.zUncorMinimizer(:,:,i), obj.GMM_MSC_zUncor(i)] = Select_IVUncorrelatedOnes(obj,lambda);
            end

            % getSelectedIdentities
            % For IV Uncorrelated Part
            [~,index] = min(obj.GMM_MSC_zUncor);
            obj.SelectedLambdazUncor = obj.LambdaRange(index);

            %obj.SelectedzUncor  = transpose(double(all(obj.zUncorMinimizer(:,:,index)==0)));
            obj.SelectedzUncor  = double(obj.zUncorMinimizer(:,:,index)==0);

            obj.Num_Selected_Candidate_Control = sum(obj.SelectedzUncor);
        end

        %% Auxiliary Public Functions 
        function value = isAllTrue_z_UncorSelected(obj)
            if all(obj.SelectedzUncor == obj.True_z_Uncor)
                value = 1;
            else
                value = 0;
            end
        end

        function  SelectedControlData = getSelectedControlData(obj)
            SelectedControlData = obj.candidateData(:,logical(obj.SelectedzUncor));
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
    end

    %% Methods For Simulations
    methods(Access = public)
        function FinalizeSelection(obj)
            obj.zData = [];
            obj.candidateData = [];
        end
    end
end
