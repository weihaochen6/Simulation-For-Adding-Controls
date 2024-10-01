classdef iidErrorOnlyAddControlsDataGenerator<handle
    %% INPUTS Properties
    properties(Access = public)
        %% Basic Data Dimensions
        T (1,1){mustBePositive,mustBeInteger} = 100

        %% Parameters setting
        %% coefficients
        % theta is the coefficient of x in y = x * Theta0 + u
        Theta0 = 0
        % Gamma is the coefficient of z in x = w' * Gamma + v
        Gamma_1 = 0.1
        Gamma_2 = 0.1 
        
        Gamma_3 = 0

        Alpha_z = 0.2
        Alpha_f = 0.2

        Pi_z = 0.2
        Pi_u = 0.2

        c_w2 = 0.2

        Num_Redundant = 1
    end

    %% OUTPUTS Properties
    properties(GetAccess = public, SetAccess = private)
        YData = []
        XData = []
        InitialValidIVData = []
        ValidControlData = []
        InvalidControlData = []
        RedundantControlData = []

        %% The True Identities Vectors
        True_z_Uncor
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
        function obj = iidErrorOnlyAddControlsDataGenerator(inputStruct)
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
            w1 = randn(obj.T,1);
           
            z = randn(obj.T,1);
            f_x = randn(obj.T,1);
            v =  randn(obj.T,1);
            x = z*obj.Alpha_z + obj.Alpha_f*f_x + v;

            w2_star = randn(obj.T,1);
            w2 = w2_star + obj.c_w2*f_x;
            w3 = randn(obj.T,1);

            e = randn(obj.T,1);
            u = w1*obj.Gamma_1 + w2*obj.Gamma_2 + w3*obj.Gamma_3 + e;

            y = obj.Theta0*x + u;

            w4 = obj.Pi_z * z + obj.Pi_u * u;
            w5 = randn(obj.T,obj.Num_Redundant);

            obj.YData = y;
            obj.XData = x;
            obj.InitialValidIVData = z;
            obj.ValidControlData = [w1 w2 w3];
            obj.InvalidControlData = w4;
            obj.RedundantControlData = w5;

                        %% Setting The True Identities
            obj.True_z_Uncor = [1;1;1;0;ones(obj.Num_Redundant,1)];
        end

        function CandidateData = getCandidateData(obj)
             CandidateData = [obj.ValidControlData, ...
                             obj.InvalidControlData, ...
                             obj.RedundantControlData];
        end
    end
end