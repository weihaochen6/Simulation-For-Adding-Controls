%% See Hall, 2005 Page 82 Newey and West's Method of Bandwidth Selection
%% and Newey and West, 1994
classdef NeweyWestHACEstimator <handle

    %% INPUTS Properties
    properties(Access = public)
        % T x r matrix of sample moment conditions
        SampleMoments 
        FirstIVConstant = false
    end

    %% OUTPUTS Properties
    properties(Access = public)
        LongRunVarianceEstimate = []
        OptimalBandWidth = []
    end

    %% Basic Properties
    properties(Access = public)
        T
        NumOfMoments
        WeightsOfMoments

        AHat
        HtResidual
    end


    %% Constructor Methods
    methods(Access = public)
        function obj = NeweyWestHACEstimator(SampleMoments)
            if nargin == 0
                return;
            end

            if nargin>=1 && ~isempty(SampleMoments)
                [obj.T, obj.NumOfMoments] = size(SampleMoments);
                obj.SampleMoments = SampleMoments;
                obj.WeightsOfMoments = ones(obj.NumOfMoments, 1);
            end
        end
    end

    %% Public Set Methods
    methods
        function set.FirstIVConstant(obj,value)
            if SetFirstTermConstant(value)
                obj.FirstIVConstant = value;
            end
        end

        function set.SampleMoments(obj,value)
            if SetSampleMoments(value)
                obj.SampleMoments = value;
            end
        end
    end

    %% Private Set Methods
    methods(Access = private)
        function isValid = SetFirstTermConstant(obj,value)
            if ~islogical(value)
                error('FirstIVConstant must be a logical value.');
            end

            isValid = true;
            obj.WeightsOfMoments(1) = 0;
        end

        function isValid = SetSampleMoments(obj,SampleMoments)
            if ~ismatrix(value)
                error('Sample Moments must be a matrix.');
            end

            isValid = true;
            [obj.T, obj.NumOfMoments] = size(SampleMoments);
            obj.WeightsOfMoments = ones(obj.NumOfMoments, 1);
            obj.WeightsOfMoments(1) = 1 - double(obj.FirstIVConstant);
        end
    end
    %% Main Public Methods
    methods(Access = public)
        function LongRunVarianceEstimator = getLongRunVarianceEstimate(obj)

            getOptimalBandWidth(obj);

            getUnscaledEstimator(obj);

            getScaledEstimator(obj);

            LongRunVarianceEstimator = obj.LongRunVarianceEstimate;
        end
    end

    %% Key Caclulation Steps
    methods(Access = private)
        %% See (2.2) in Newey and West, 1994
        function getOptimalBandWidth(obj)
            w = obj.WeightsOfMoments;
            Ht = obj.SampleMoments(2:end, :);% (T-1)xNumOfMoments
            HtLag= Ht(1:end-1, :);  %(T-1)xNumOfMoments

            SumHtHtLag = Ht'*HtLag;
            SumHtLagHtLag = HtLag' * HtLag;

            obj.AHat = SumHtHtLag/SumHtLagHtLag;

            obj.HtResidual = Ht - HtLag*obj.AHat;
            weighted_residual = (w')*(obj.HtResidual'); % 1 x (T-1)

            
            s0 = 0;
            s1 = 0;
            sigma0 = (1/(obj.T-1))*(weighted_residual')*weighted_residual;
            n = round(4*(obj.T/100)^(2/9));
            s0 = s0 + sigma0;
            for j = 1:n
                weighted_residual_t = (w')*(obj.HtResidual(j+1:end,:)');
                weighted_residual_tMinusj = (w')*(obj.HtResidual(1:end-j,:)');
                Sigmaj =(weighted_residual_t*(weighted_residual_tMinusj'))/(obj.T-1);
                s1 = s1 + 2*j*Sigmaj;
                s0 = s0 + 2*Sigmaj;
            end

            gamma = 1.1447*((s1/s0)^(2/3));
            obj.OptimalBandWidth = round(gamma*(obj.T^(1/3)));
        end

        function getUnscaledEstimator(obj)
            if obj.OptimalBandWidth~=0
                m = obj.OptimalBandWidth;
                Omega0 = (1/(obj.T-1))*(obj.HtResidual')*(obj.HtResidual);
                SBarStar = Omega0;
                for j = 1:m
                    Omegaj = ((obj.HtResidual(j+1:end,:)')*obj.HtResidual(1:end-j,:))/(obj.T-1);
                    % Use Bartlett kernel
                    kj= 1-j/(m+1);
                    SBarStar=SBarStar+kj*(Omegaj+Omegaj');
                end
                weight = eye(obj.NumOfMoments)/(eye(obj.NumOfMoments) - obj.AHat);
                obj.LongRunVarianceEstimate = weight*SBarStar*(weight');
            end
        end
    end
end