classdef Rothfelder_Boldea2022<handle
    %% INPUTS Properties
    properties(Access = public)
        T = 100

        Delta_x = 0
        Gamma0 = 2.25
        Delta_Pi = 0.5
        Rho0 = 1.75

        Variance_u = 1
        Variance_e = 1
        CovarianceOfeAndu = 0

        ErrorType = 'homoskedastic'
    end
    %% OUTPUTS Properties
    properties(GetAccess = public, SetAccess = private)
        yData = []
        xData = []
    end

    %% String Properties For Estimation Specification
    properties(Constant,Access = private)
        HomoskedasticStr = 'homoskedastic'
        ConditionalHeteroskedastic= 'conditional heteroskedastic'
    end

    %% Constructor Methods
    methods(Access = public)
        function obj = Rothfelder_Boldea2022(inputStruct)
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
        function obj = generateData(obj)
            % generate  zt
            z_t = 1 + randn(obj.T,1);
            q_t = z_t + 1;
            mu = [0 0];
            Sigma = eye(2);
            Sigma(1,2) = obj.CovarianceOfeAndu;
            Sigma(2,1) = obj.CovarianceOfeAndu;
            Data = mvnrnd(mu,Sigma, obj.T);

            e_t = Data(:,1);
            u_t = Data(:,2);

            if strcmp(obj.ErrorType,obj.HomoskedasticStr)
                e_t = randn(obj.T,1);
                epsilon_t = e_t;

            elseif strcmp(obj.ErrorType,obj.ConditionalHeteroskedastic)

                epsilon_t = e_t*z_t/sqrt(2);
            else
                error('Wrong Input For Error Type.');

            end
            x_t = 1 + z_t + obj.Delta_Pi*z_t*(q_t>obj.Rho0) + u_t;
            y_t = 1 + x_t + obj.Delta_x*(q_t>obj.Gamma0) + epsilon_t;

            obj.xData = x_t;
            obj.yData = y_t;
        end
    end
end

