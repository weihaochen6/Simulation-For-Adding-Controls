classdef iidErrorCaseDataGenerator<handle
    %% INPUTS Properties
    properties(Access = public)
        %% Basic Data Dimensions
        T (1,1){mustBePositive,mustBeInteger} = 100

        %% Parameters setting
        %% coefficients
        % theta is the coefficient of x in y = x * Theta0 + u
        Theta0 = 0
        % Gamma is the coefficient of z in x = z' * Gamma + v
        Gamma_S = 0.1
        Gamma_A0 = 0.1
        Gamma_A1 = 0.1 
        Gamma_W = 0.2

        c_zA1 = 0.1
        c0 = 0.2
        cbar = 2.4
        c0_B2 = 0.2 
        cbar_B2 = 2.4
        % The Total Number of Candidate Set of IVs and Controls
        K = 7 
    end

    %% OUTPUTS Properties
    properties(GetAccess = public, SetAccess = private)
        YData = []
        XData = []
        InitialValidIVData = []
        RelevantIV_InvalidControlData = []
        RelevantIV_ValidControlData = []
        RedundantIV_RedundantControlData = []
        InvalidIV_RedundantControlData = []
        InvalidIV_InvalidControlData = []
        InvalidIV_RelevantControlData = []

        %% The True Identities Vectors
        True_zS_Uncor
        True_u_Uncor

        True_Add_As_IV
        True_Add_As_Contol
    end

    %% Basic Properties
    properties(GetAccess = public, SetAccess = private)
        Num_RedunIV_RedunControl
        Num_InvalidIV_ValidControl
        Num_InvalidIV_InvalidControl
    end

    %% Constructor Methods
    methods(Access = public)
        function obj = iidErrorCaseDataGenerator(inputStruct)
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

    %% Set Methods
    methods
        function  set.K(obj, value)
            if mod(value - 3, 4) == 0 
                obj.K = value;
            else
                error('The value of K must be a multiple of 4 plus 3.');
            end
        end
    end

    %% Main Public Methods
    methods(Access = public)
        function obj = generateData(obj)
            %% Create The Total Covariance and Variance Matrix
            Sigma_AS = eye(3);
            Sigma_AS(1,3) = 0.2;
            Sigma_AS(3,1) = 0.2;
            Sigma_B = eye(obj.K-3);
            Sigma_w = 1;
            Sigma_ev = [0.5 0.6; 0.6 1];
            Sigma = blkdiag(Sigma_AS, Sigma_B, Sigma_w, Sigma_ev);
            mu = zeros(obj.K + 3,1);
            Data = mvnrnd(mu,Sigma, obj.T);

            z_S = Data(:,1);
            z_A0 = Data(:,2);
            z_A1 = Data(:,3);

            obj.Num_RedunIV_RedunControl = (obj.K-3)/2;
            index = 3;
            z_B0 = Data(:,index+1:index + obj.Num_RedunIV_RedunControl);
            index = index + obj.Num_RedunIV_RedunControl;
            
            obj.Num_InvalidIV_ValidControl = (obj.K-3)/4;
            z_B1Star = Data(:,index+1:index+obj.Num_InvalidIV_ValidControl);
            index = index+obj.Num_InvalidIV_ValidControl;

            obj.Num_InvalidIV_InvalidControl = (obj.K-3)/4;
            z_B2Star = Data(:,index+1:index + obj.Num_InvalidIV_InvalidControl);
            index = index+obj.Num_InvalidIV_InvalidControl;
            
            w = Data(:,index+1);
            index = index + 1;

            e = Data(:,index+1);
            index= index + 1;
            v = Data(:,index+1);

            %% Obtain u, x, y zB1, zB2
            u = obj.Gamma_W*w + e;
            z_B1 = zeros(size(z_B1Star));
            z_B2 = zeros(size(z_B2Star));
            for l = 1:obj.Num_InvalidIV_ValidControl
                c_ell = EndogeneityStrength(l,obj.c0,obj.cbar,obj.Num_InvalidIV_ValidControl);
                z_B1(:,l)= z_B1Star(:,l) + obj.c_zA1*z_A1 + c_ell*u;
            end

            for l = 1:obj.Num_InvalidIV_InvalidControl
                c_ell = EndogeneityStrength(l,obj.c0,obj.cbar,obj.Num_InvalidIV_InvalidControl);
                c_ell_B2 = EndogeneityStrength(l,obj.c0_B2,obj.cbar_B2,obj.Num_InvalidIV_InvalidControl);
                z_B2(:,l)= z_B2Star(:,l) + c_ell*u + c_ell_B2*z_S;
            end

            x = z_S*obj.Gamma_S + z_A0*obj.Gamma_A0 + z_A1*obj.Gamma_A1 + v;
            y = x*obj.Theta0 + u;

            %% Store The Data
            obj.YData = y;
            obj.XData = x;
            obj.InitialValidIVData = z_S;
            obj.RelevantIV_InvalidControlData = z_A0;
            obj.RelevantIV_ValidControlData = z_A1;
            obj.RedundantIV_RedundantControlData = z_B0;
            obj.InvalidIV_RedundantControlData = z_B1;
            obj.InvalidIV_InvalidControlData = z_B2;
            obj.InvalidIV_RelevantControlData = w;

            %% Setting The True Identities
            obj.True_zS_Uncor = [1;...
                                 0;...
                                 ones(obj.Num_RedunIV_RedunControl,1); ...
                                 ones(obj.Num_InvalidIV_ValidControl,1);...
                                 zeros(obj.Num_InvalidIV_InvalidControl,1);...
                                 1];


            obj.True_u_Uncor = [1;...
                                1;...
                                ones(obj.Num_RedunIV_RedunControl,1);...
                                zeros(obj.Num_InvalidIV_ValidControl,1);...
                                zeros(obj.Num_InvalidIV_InvalidControl,1);...
                                0];

            obj.True_Add_As_IV = [1;...
                                  1;...
                                  zeros(obj.K-2,1)];
            
            obj.True_Add_As_Contol = [zeros(obj.K-1,1);1];
        end

        function CandidateData = getCandidateData(obj)
             CandidateData = [obj.RelevantIV_InvalidControlData, ...
                             obj.RelevantIV_ValidControlData, ...
                             obj.RedundantIV_RedundantControlData, ...
                             obj.InvalidIV_RedundantControlData, ...
                             obj.InvalidIV_InvalidControlData, ...
                             obj.InvalidIV_RelevantControlData];
        end
    end
end

function c_ell = EndogeneityStrength(ell,c0,cbar,Num)
c_ell = c0 + (ell-1)*(cbar-c0)/Num;
end