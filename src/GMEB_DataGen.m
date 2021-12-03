%   File: GMEB_DataGen.m
%   Copyright (c) <2019> <University of Mons>
%   COLORAMAP Group, Univ. of Mons, Belgium
%   https://sites.google.com/site/nicolasgillis/projects/overview
%
%   Permission is hereby granted, free of charge, to any person
%   obtaining a copy of this software and associated documentation
%   files (the "Software"), to deal in the Software without restriction,
%   including without limitation the rights to use, copy, modify and
%   merge the Software, subject to the following conditions:
%
%   1.) The Software is used for non-commercial research and
%       education purposes.
%
%   2.) The above copyright notice and this permission notice shall be
%       included in all copies or substantial portions of the Software.
%
%   3.) Publication, Distribution, Sublicensing, and/or Selling of
%       copies or parts of the Software requires special agreements
%       with the University of Mons and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Mons is granted 
%       the non-exclusive right to publish modifications or contributions 
%       in future versions of the Software free of charge.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%   OTHER DEALINGS IN THE SOFTWARE.
%
%   Persons using the Software are encouraged to notify the
%   COLORAMAP Group, Univ. of Mons, Belgium about bugs. Please reference 
%   the Software in your publications if it was used for them.
% ------------------------------------------------------------------------
% SYNTAX:  
% [ data, param ] = GMEB_DataGen( param );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% Generates points sampled uniformly from a ball on a Grassmann manifold.
%
% Includes options to:
% a.) sample uniformly from the boundary of the ball
% b.) sample from the boundary of a second, smaller baller nested within 
%       the larger one
% c.) sample non-uniformly from the boundary of the ball (with a more 
%       densely sampled arc), 
% d.) sample from the interior of the aforementioned balls. 
% 
% The samples form a collection of data whose for which the center of the 
% minimum enclosing ball may not equal to the center of mass.
%
% Random additional dimensions and white Gaussian noise are then added to
% create data on a collection of Grassmann manifolds for which the optimal
% minimal enclosing ball dimension is known.
%
% ------------------------------------------------------------------------
% INPUTS:
% param     MATLAB structure that contains the parameters for data
%           generation. Minimally includes the field 'default' with value
%           [true].
%
%   Fields and [default values]:
%   'default'           boolean - required
%                       [true]
%   'ambient_dims'      scalar - dimension of the ambient space
%                       [100]
%   'num_samples_bb'    scalar - # of samples from big ball
%                       [10]
%   'num_samples_sb'    scalar - # of samples from small ball
%                       [10]
%   'common_dims_bb'    scalar - subspace dimensions shared by all big ball 
%                       points 
%                       [3]
%   'common_dims_sb'    scalar - subspace dimensions shared by all small 
%                       ball points
%                       [3]
%   'max_dims'          scalar - maximum dimension of any subspace
%                       [10]
%   'big_rad'           scalar - radius of big ball
%                       [.5]
%   'center_sep'        scalar - distance between center of big ball and 
%                       center of small ball
%                       [(.75)*big_rad]
%   'small_rad'         scalar - radius of small ball
%                       [(.25)*center_sep]
%   'noise_var'         scalar - variance of the white Gaussian noise
%                       [0.001]
%   'use_subsample'     boolean - flag to generate additional samples for
%                       the big ball and subsample to increase uniformity
%                       [false]
%   'subsample_percent' scalar - percent of total sample to keep as the
%                       subsample (subsample will always have
%                       num_samples_bb points)
%                       [.5]
%   'eq_variance'       boolean - if true, added subspace dimensions will
%                       include signal of equal variance to correlated 
%                       signals.  otherwise, the extra dimensions are
%                       noise. If true, noise_var = 0; means that there is
%                       no additive gaussian noise in noisy samples but
%                       there is extra dimensions with per component
%                       variance equal to signal dimensions.
%                       [true]
%   'dims_bb'           num_samples_bb x 1 double - dimensions of big ball 
%                       points after additional orthogonal dimensions
%                       have been added
%                       [randi([common_dims_bb,max_dims],[num_samples_bb,1])]
%   'dims_sb'           num_samples_sb x 1 double - dimensions of small
%                       ball points after additional orthogonal dimensions
%                       have been added
%                       [randi([common_dims_sb,max_dims],[num_samples_sb,1])]
%   'nonuni'            boolean - if true, extra points will be generated
%                       nonuniformly in the neighborhood of one point on 
%                       the boundary of the large ball.  otherwise, the 
%                       points are approximately uniform
%                       [false]
%   'nonuni_amt'        integer - selects the number of extra points that 
%                       come from the nonuniform region
%                       [5]
%   'd_type'            string - specifies the distance type to be used.
%                       options are 'DPF' for the squared projection
%                       Frobenious norm or 'DGEO' for the geodesic
%                       distance based on arc-length.
%                       ['DPF']
%   'int_pts'           boolean - if true, points will be generated
%                       uniformly within the interior of both balls in the
%                       same ratio as the boundary poitns.  otherwise, the 
%                       points are all on the boundary
%                       [false]
%   'int_pts_ratio'     scalar - selects the number of points on the
%                       interior of each ball to be ceil(int_pts_ratio * 
%                       num_samples_bb) (or num_samples_sb respectively).\
%                       [1]
%   'dims_bb_int'       ceil(int_pts_ratio * num_samples_bb) x 1 double - 
%                       dimensions of big ball interior points after 
%                       additional orthogonal dimensions have been added
%                       [randi([common_dims_bb,max_dims],
%                       [ceil(int_pts_ratio * num_samples_bb),1])]
%   'dims_sb_int'       ceil(int_pts_ratio * num_samples_sb) x 1 double - 
%                       dimensions of small ball interior points after 
%                       additional orthogonal dimensions have been added
%                       [randi([common_dims_sb,max_dims],
%                       [ceil(int_pts_ratio * num_samples_sb),1])]
%
% ------------------------------------------------------------------------
% OUTPUTS: 
% data      MATLAB structure with fields for data points computed GMEB
%           solution, if requested.
%   Fields:     
%   center:             ambient_dims x common_dims_bb orthonormal matrix         
%   clean_samples_bb:   num_samples_bb x 1 cell array
%   clean_samples_sb:   num_samples_sb x 1 cell array
%   clean_samples_all:  (num_samples_bb + num_samples_sb) x 1 cell array
%   noisy_samples_bb:   num_samples_bb x 1 cell array
%   noisy_samples_sb:   num_samples_sb x 1 cell array
%   noisy_samples_all:  (num_samples_bb + num_samples_sb) x 1 cell array
%
%   These fields are only included if 'interior_pts' = 'TRUE'
%   clean_samples_bb_int:   ceil(int_pts_ratio*num_samples_bb) x 1 cell
%                           array
%   clean_samples_sb_int:   ceil(int_pts_ratio*num_samples_sb) x 1 cell 
%                           array
%   clean_samples_all_int:  (ceil(int_pts_ratio*num_samples_bb) + ceil(
%                           int_pts_ratio*num_samples_sb)) x 1 cell array
%   noisy_samples_bb_int:   ceil(int_pts_ratio*num_samples_bb) x 1 cell
%                           array
%   noisy_samples_sb_int:   ceil(int_pts_ratio*num_samples_sb) x 1 cell 
%                           array
%   noisy_samples_all_int:  (ceil(int_pts_ratio*num_samples_bb) + ceil(
%                           int_pts_ratio*num_samples_sb)) x 1 cell array
%   clean_samples_both:     (num_samples_bb + num_samples_sb + 
%                           ceil(int_pts_ratio*num_samples_bb) + 
%                           ceil(int_pts_ratio*num_samples_sb)) x 1 cell 
%                           array
%   noisy_samples_both:     (num_samples_bb + num_samples_sb + 
%                           ceil(int_pts_ratio*num_samples_bb) + 
%                           ceil(int_pts_ratio*num_samples_sb)) x 1 cell 
%                           array
%
% param     MATLAB structure with the parameters used for data generation.
%           Only differs from input structure if fields were left blank.
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [1] SetParam.m
%       Included within.
%       Sets default parameters.
%       Written by T. Marrinan; see file for documentation.
%
% [2] GrGeodesic.m
%       Included within.
%       Returns a point along a Grassmannian endpoint geodesic.
%       Written by T. Marrinan; see file for documentation.
%
% [3] GrDist.m
%       Included within.
%       Computes a variety of Grassmannian metrics.
%       Written by T. Marrinan; see file for documentation.
%
% ------------------------------------------------------------------------
% DETAILS:
% The purpose of this data generation function is to create samples for 
% which the ground truth solution for the center of the minimum enclosing 
% ball is known on Gr(k,n), the Grassmann manifold of k-planes in R^n.
%
% To have that, I propose sampling points uniformly distributed on a unit 
% ball on the Grassmannian Gr(k,n) for a chosen k and then completing these
% points to bases on Gr(p_i,n) for p_i > k.
%
% A second, smaller ball is then centered completely within the previous
% one and points are sampled from within the smaller ball.
%
% The union of those two point sets creates a data set that has a different
% karcher mean and gmeb center, both of which can be found empirically.
%
% The method for creating this data is as follows:
%
% (1) Sample P < M points uniformly on Gr(k,n) using the method of [2].
%     a. Generate n x n matrix, S, with Gaussian entries.
%     b. Compute S = QR with any method.
%     c. Compute T = Q*diag(sign(diag(R))) as suggested by [3].
%     d. Resample the initial random draw using the method of [4].
%        i. Let D(x) denote the shortest distance from a data point x to 
%           the closest center we have already chosen. This means (for M
%           points), M-i new distances need to be computed at each step
%       ii. Choose an initial center c_1 uniformly at random from X.
%      iii. Choose the next center ci, selecting c_i = x' ? X with 
%           probability D(x')^2/ Sum_{x?X}( D(x)^2 ).
%       iv. Repeat Step iii. until we have chosen a total of P centers.
%
% (2) Compute the extrinsic mean, mu, of the P points.
%
% (3) Estimate a uniform sample on the unit ball around [mu]. That is,
% take a geodesic step from mu in the direction of each of the P points,
% c_i such that d(mu,c_i) = a for some fixed a>0. The step depends on the
% distance measure, and is currently possible for the projection
% Frobenius norm (squared) and for the geodesic distance based on arc
% -length.
%
% (4) Recompute the extrinsic mean of the points on the sphere, if they
% are indeed unifom, the newly compute mean should be the same as the
% previous one, right?
%
% Option b.)
% (5) Compute a point, x*, some distance b>0 between the center and a 
% random point on the boundary.
%
% (6) Generate another uniform sampling of points and map them to an
% epsilon ball around x*.
%
% Option c.)
% (7) Randomly select two points on the boundary of the large ball,
% [endPt1] and [endPt2].
% 
% (8) Sample points uniformly at random from the minimal geodesic
% connecting the endpoints.
%
% (9) Move the newly sampled points to the boundary of the large sphere by
% taking a geodesic step from the center, mu, in the direction of each of
% these new points.  The result should be points with distance a from the
% center, that are only present within one arc of the boundary.
%
% Option d.)
% (10) Points can then be added to the interior of either ball by taking
% steps of length less than a or b from the centers of the respective
% balls
% 
% Final steps
% (11) Including all the new points and averaging all the points should move
% the Karcher mean, flag mean, and geometric mean, but not the GMEB center.
% This should be a proof of concept that things are being computed
% correctly.
%
% (12) Complete bases of the existing points to higher dimensions and run
% the optimal dimension stuff.
%
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% [2] 	Liu, Shusen, P.T. Bremer, J. J. Jayaraman, Bei Wang, Brian Summa, 
% 	and Valerio Pascucci. "The Grassmannian Atlas: A General Framework for 
% 	Exploring Linear Projections of High Dimensional Data." In Computer 
% 	Graphics Forum, vol. 35, no. 3, pp. 1-10. 2016.
%
% [3] 	Mezzadri, Francesco. "How to generate random matrices from the 
% 	classical compact groups." arXiv preprint math-ph/0609050 (2006).
%
% [4] 	Arthur, David, and Sergei Vassilvitskii. "k-means++: The advantages 
% 	of careful seeding." In Proceedings of the eighteenth annual ACM-SIAM 
% 	symposium on Discrete algorithms, pp. 1027-1035. Society for Industrial 
% 	and Applied Mathematics, 2007.
%
% [5] 	Renard, Emilie, Kyle A. Gallivan, and P-A. Absil. "A Grassmannian 
% 	Minimum Enclosing Ball Approach for Common Subspace Extraction." In 
% 	International Conference on Latent Variable Analysis and Signal 
% 	Separation, pp. 69-78. Springer, Cham, 2018.
%
% [6] 	Gallivan, Kyle A., Anuj Srivastava, Xiuwen Liu, and Paul Van Dooren. 
% 	"Efficient algorithms for inferences on grassmann manifolds." In IEEE 
% 	Workshop on Statistical Signal Processing, 2003, pp. 315-318. IEEE, 2003.
% 
% ------------------------------------------------------------------------
% CREATED:      20/05/2019 by Tim Marrinan
%
% LAST EDITED:  11/03/2020 by Tim Marrinan
%
% NOTES: 
% 22/05/2019:
% Updated the GrGeodesic function to use the appropriate estimate when the
% principal angle is close to zero. See [6] for details. (TM) 
%
% 03/10/2019:
% Updated to generate interior points in both balls as well with the same
% ratio as the boundary points. (TM)
%
% 16/10/2019:
% Updated so that random additional dimensions are orthogonal to center &
% current dimensions if possible. Only affects points on the boundary of
% the large ball. Makes the center more accurate. (TM)
%
% 11/03/2020:
% A unit ball varies depending on the choice of distance function, so that
% needs to be a parameter as well. (not included yet) (TM)
% ------------------------------------------------------------------------
function [ data, param ] = GMEB_DataGen( param )

    %% Set parameters
    param = SetParam( param ); % deprecated default parameter function
    % More condensed variable names because I am lazy
    n                   = param.ambient_dims;
    nBB                 = param.num_samples_bb;
    nSB                 = param.num_samples_sb;
    k                   = param.common_dims_bb;
    small_k             = param.common_dims_sb;
    kmax                = param.max_dims;
    dims_bb             = param.dims_bb;
    dims_sb             = param.dims_sb;
    big_rad             = param.big_rad;
    center_sep          = param.center_sep;
    small_rad           = param.small_rad;
    sigmaN              = param.noise_var;
    useSubsample        = param.use_subsample;
    subsamplePercent    = 1/param.subsample_percent;
    eqVariance          = param.eq_variance;
    nonuni              = param.nonuni;
    nonuni_amt          = param.nonuni_amt;
    int_pts             = param.int_pts;
    d_type              = param.d_type; % 'DPF'; % 'DGEO';
    if int_pts
        int_pts_ratio   = param.int_pts_ratio;
        dims_bb_int     = param.dims_bb_int;
        dims_sb_int     = param.dims_sb_int;
    end
    
    
    
    % If there are no common dimensions in the data, everything else is
    % irrelevant.  We just need to compute a uniform sample on a colletion
    % of Grassmannians.
    if param.common_dims_bb == 0
        % If there are no common dimensions, we have to generate the
        % data a little differently.  This should probably get added to
        % the data generation function.
        n = param.ambient_dims;
        kmax = param.max_dims;
        M = param.num_samples_bb;
        S = randn(n,kmax,M);
        T = zeros(size(S));
        Q = zeros(size(S));
        R = zeros(kmax,kmax);
        data.noisy_samples_all = cell(M,1);
        for i = 1 : M
            [Q(:,:,i), R(:,:,i)] = qr(S(:,:,i),0);
            T(:,:,i) = Q(:,:,i)*diag(sign(diag(R(:,:,i))));
            param.dims_bb = randi([kmax-2,param.max_dims],...
                [param.num_samples_bb,1]);
            param.dims_sb = randi([kmax-2,param.max_dims],...
                [param.num_samples_sb,1]);
            data.noisy_samples_all{i} = T(:,1:param.dims_bb(i),i);
        end
        return;
    end


    %% Generate uniformly distributed orthonormal matrices, via [2]
    if useSubsample
        M = ceil(subsamplePercent*nBB);
    else
        M = nBB;
    end

    S = randn(n,k,M);
    T = zeros(size(S));
    Q = zeros(size(S));
    R = zeros(k,k);
    for i = 1 : M
        [Q(:,:,i), R(:,:,i)] = qr(S(:,:,i),0);
        T(:,:,i) = Q(:,:,i)*diag(sign(diag(R(:,:,i))));
    end
    X = T(:,1:k,:);


    %% Compute extrinsic mean, MU, of all the data.
    % It is perhaps more appropriate to compute the Karcher mean here, but
    % practically it is unnecessary because the points need to be moved to
    % a fixed distance away from the center anyway.  This 'MU' is just an
    % so that not all points end up on on portion of the ball.
    uniform_samples = cell(M,1);
    for i = 1 : M
        uniform_samples{i} = X(:,:,i);
    end
    [flag,~,~] = svd([uniform_samples{:}],'econ');
    MU = flag(:,1:k);
    data.center = MU;

    %% Project the points onto the boundary of a ball around MU
    unitBall_samples = cell(M,1);
    for i = 1 : M
        old_dist = GrDist(MU,uniform_samples{i},k,d_type)^2;
        unitBall_samples{i} = GrGeodesic(MU,uniform_samples{i},...
            k,old_dist,big_rad,d_type); % Defined below.
    end
    

    % ------------------------------------------------------------------------
    % NOTE:
    % Perhaps I should run the resampling from [3] after mapping the 
    % points to the circle rather than before? There is no guarrantee that the
    % regularity from the initial sample is preserved after shrinking the
    % distance is there?  Especially when all the points are outside of the
    % convex radius. This is the current method. 
    % ------------------------------------------------------------------------


    %% Resample points to improve regularity via [1] using the method of [3]
    if useSubsample
        full_ids = 1:M;
        c_ids = false(M,1);
        new_center = randi(M);
        c_ids(new_center) = true;
        C = unitBall_samples{new_center};
        min_dists = 100+zeros(M,1);
        min_dists = min_dists(~c_ids);
        for i = 1 : nBB-1
            remaining_ids = full_ids(~c_ids);
            for j = 1 : M-sum(c_ids)
                temp_dist = GrDist(C,unitBall_samples{...
                    remaining_ids(j)},k,d_type); % Defined below
                min_dists(j) = min(min_dists(j),temp_dist^2);
            end
            bins = cumsum(min_dists./sum(min_dists));
            draw = rand;
            sample = find(cumsum(bins)>draw,1);
            new_center = remaining_ids(sample);
            c_ids(new_center) = true;
            C = unitBall_samples{new_center};
            min_dists(sample) = [];
        end
        unitBall_samples = unitBall_samples(c_ids);
    end
      
    %% Average and center again, because the subsampling may change things
    [flag,~,~] = svd([unitBall_samples{:}],'econ');
    MU = flag(:,1:k);
    data.center = MU;
    for i = 1 : nBB
        old_dist = GrDist(MU,unitBall_samples{i},k,d_type)^2;
        unitBall_samples{i} = GrGeodesic(MU,uniform_samples{i},...
            k,old_dist,big_rad,d_type); % Defined below.
    end
    
    %% Generate extra samples on the big ball in the neighborhood of one pt
    if nonuni
        ids = randperm(nBB,2);
        nonuni_samples = cell(nonuni_amt,1);
        % -----------------------------------------------------------------
        % NOTE:
        % 29.11.2019 I am going to make it possible for the high density 
        % area have a higher common dimenison
        % -----------------------------------------------------------------
        endpt_1 = unitBall_samples{ids(1)};
        endpt_2 = unitBall_samples{ids(2)};
        between_dist = GrDist(endpt_1,endpt_2,k,'DGEO');
        if small_k > k
            extra_dims  = randn(n,small_k-k);
            [extra_dims, ~] = qr(extra_dims- MU*(MU'*extra_dims),0);
        end   
        for i = 1 : nonuni_amt
            % Put a point a random percentage of the distance between the
            % two random arc endpoints.
            % (the distance here should be geodesic distance because we
            % just want a more or less uniform distribution on the arc)
            nonuni_samples{i} = GrGeodesic(endpt_1,endpt_2,k,...
                between_dist,rand*between_dist,'DGEO'); % Defined below.
            old_dist = GrDist(MU,nonuni_samples{i},k,d_type)^2; % This should only be squared if that is specified somewhere
            % Readjust the distance between the new point and the center so
            % that the new point is actually on the boundary of the ball.
            nonuni_samples{i} = GrGeodesic(MU,nonuni_samples{i},...
                k,old_dist,big_rad,d_type); % Defined below.
            % -------------------------------------------------------------
            % NOTE:
            % 29.11.2019 I am going to make it possible for the high 
            % density area have a higher common dimenison
            % -------------------------------------------------------------
            if small_k > k
                nonuni_samples{i} = [nonuni_samples{i},extra_dims];
                [nonuni_samples{i},~,~] = svd(nonuni_samples{i},'econ');
                nonuni_samples{i} = nonuni_samples{i}(:,1:small_k);
            end
        end
        unitBall_samples = [unitBall_samples;nonuni_samples];
        param.num_samples_bb = nBB+nonuni_amt;
        nBB = nBB+nonuni_amt;
        % -----------------------------------------------------------------
        % NOTE:
        % 29.11.2019 I am going to make it possible for the high density 
        % area have a higher common dimenison
        % -----------------------------------------------------------------
        param.dims_bb = [param.dims_bb;...
            randi([param.common_dims_sb,param.max_dims],...
            [param.nonuni_amt,1])];
        dims_bb = param.dims_bb;
    end

    data.clean_samples_bb = unitBall_samples;

    %% Generate points from a second uniform distribution
    % It's not so important that these have the same regularity is the 
    % initial draw so skip the resampling.
    S2 = randn(n,small_k,nSB);
    T2 = zeros(size(S2));
    Q2 = zeros(size(S2));
    R2 = zeros(small_k,small_k);
    for i = 1 : nSB
        [Q2(:,:,i), R2(:,:,i)] = qr(S2(:,:,i),0);
        T2(:,:,i) = Q2(:,:,i)*diag(sign(diag(R2(:,:,i))));
    end
    X2 = T2(:,1:small_k,:);


    %% Project these new points into a neighborhood within the unit ball
    % Find a point between the the center and the boundary of the ball.
    % Project into an epsilon neighborhood of that point
    small_center_id = randi(nBB);
    small_center_direction = unitBall_samples{small_center_id};
    data.small_center_id = small_center_id;
    small_center = GrGeodesic(MU,small_center_direction,k,big_rad,...
        center_sep,d_type); % Defined below.
    smallBall_samples = cell(nSB,1);
    small_samples = cell(nSB,1);
    if small_k > k
        clear temp_sc u s v
        temp_sc = [small_center randn(n,small_k-k)];
        [temp_sc,~,~] = svd(temp_sc,'econ');
        temp_sc = temp_sc(:,1:small_k);
    else
        temp_sc = small_center(:,1:small_k);
    end   
    for i = 1:nSB
        small_samples{i} = X2(:,:,i);     
        old_dist = GrDist(temp_sc,small_samples{i},small_k,d_type)^2;
        smallBall_samples{i} = GrGeodesic(temp_sc,small_samples{i},...
            small_k,old_dist,small_rad,d_type); % Defined below.
    end
    data.small_center = temp_sc;
    data.clean_samples_sb = smallBall_samples;
    data.clean_samples_all = [data.clean_samples_bb; data.clean_samples_sb];
    
    %% Add a random number of dimensions to each point
    % eqVariance = true means that the "pure noise" dimensions are unit vectors
    % eqVariance = false means that the "pure noise" dimensions are vectors
    % of length sigmaN / (max(dims_sb(i),small_k))
    if eqVariance
        S3 = randn(n,kmax,nBB);
        T3 = zeros(size(S3));
        Q3 = zeros(size(S3));
        R3 = zeros(kmax,kmax,nBB);
        sloppyBall_samples = cell(nBB,1);
        for i = 1 : nBB
            [oop,~] = qr([unitBall_samples{i}, data.center],0);
            [Q3(:,:,i), R3(:,:,i)] = qr(S3(:,:,i)- oop*(oop'*S3(:,:,i)),0); % 16/10/2019 project new dimensions away from old & center
            T3(:,:,i) = Q3(:,:,i)*diag(sign(diag(R3(:,:,i))));
            % -------------------------------------------------------------
            % NOTE:
            % 29.11.2019 updated the noise variance to be overall rather
            % than per component. More components means less noise per 
            % component
            % -------------------------------------------------------------
            temp = [unitBall_samples{i}, T3(:,1:dims_bb(i)-k,i)]...
                + sqrt(sigmaN)*(1/(sqrt(dims_bb(i))...
                *sqrt(n)))*randn(n,dims_bb(i));
            [temp,~,~] = svd(temp,'econ');
            sloppyBall_samples{i} = temp(:,1:dims_bb(i));    
        end
        clear dims S3 T3 Q3 R3

        S4 = randn(n,kmax,nSB);
        T4 = zeros(size(S4));
        Q4 = zeros(size(S4));
        R4 = zeros(kmax,kmax,nSB);
        sloppySmallBall_samples = cell(nSB,1);
        for i = 1 : nSB
            [Q4(:,:,i), R4(:,:,i)] = qr(S4(:,:,i),0);
            T4(:,:,i) = Q4(:,:,i)*diag(sign(diag(R4(:,:,i))));
            % -------------------------------------------------------------
            % NOTE:
            % 29.11.2019 updated the noise variance to be overall rather
            % than per component. More components means less noise per 
            % component
            % -------------------------------------------------------------
            temp = [smallBall_samples{i}, ...
                T4(:,1:max(0,dims_sb(i)-small_k),i)]...
                + sqrt(sigmaN)*(1/(sqrt(max(dims_sb(i),small_k))...
                *sqrt(n)))*randn(n,max(dims_sb(i),small_k));
            [temp,~,~] = svd(temp,'econ');
            sloppySmallBall_samples{i} = temp(:,1:dims_sb(i)); 
        end
    else
        T3 = randn(n,kmax,nBB);
        sloppyBall_samples = cell(nBB,1);
        for i = 1 : nBB
            % -------------------------------------------------------------
            % NOTE:
            % 29.11.2019 updated the noise variance to be overall rather
            % than per component. More components means less noise per 
            % component
            % -------------------------------------------------------------
            T3(:,:,i) = sqrt(sigmaN)*(1/( sqrt(n)...
                *sqrt(dims_bb(i))))*T3(:,:,i);
            T3(:,1:dims_bb(i),i) =  T3(:,1:dims_bb(i),i) + ...
                [unitBall_samples{i},...
                zeros(n,dims_bb(i)-size(unitBall_samples{i},2),1)];
            [temp,~,~] = svd(T3(:,1:dims_bb(i),i),'econ');
            sloppyBall_samples{i} = temp(:,1:dims_bb(i));    
        end
        clear dims S3 T3 Q3 R3

        T4 = randn(n,kmax,nSB);
        sloppySmallBall_samples = cell(nSB,1);
        for i = 1 : nSB
            % -------------------------------------------------------------
            % NOTE:
            % 29.11.2019 updated the noise variance to be overall rather
            % than per component. More components means less noise per 
            % component
            % -------------------------------------------------------------
            T4(:,:,i) = sqrt(sigmaN)*(1/(sqrt(n)*sqrt(dims_sb(i))))...
                *T4(:,:,i);
            T4(:,1:small_k,i) = T4(:,1:small_k,i) + smallBall_samples{i};
            [temp,~,~] = svd(T4(:,1:dims_sb(i),i),'econ');
            sloppySmallBall_samples{i} = temp(:,1:dims_sb(i)); 
        end
    end
    data.noisy_samples_bb = sloppyBall_samples;
    data.noisy_samples_sb = sloppySmallBall_samples;
    data.noisy_samples_all = [data.noisy_samples_bb;...
        data.noisy_samples_sb];
    
    %% Generate points in the interior of each ball
    if int_pts
        Mint = ceil(nBB*int_pts_ratio);
        Sint = randn(n,k,Mint);
        Tint = zeros(size(Sint));
        Qint = zeros(size(Sint));
        Rint = zeros(k,k);
        for i = 1 : Mint
            [Qint(:,:,i), Rint(:,:,i)] = qr(Sint(:,:,i),0);
            Tint(:,:,i) = Qint(:,:,i)*diag(sign(diag(Rint(:,:,i))));
        end
        Xint = Tint(:,1:k,:);
        uniform_samples_int = cell(Mint,1);
        for i = 1 : Mint
            uniform_samples_int{i} = Xint(:,:,i);
        end
        
        % Project these points into the interior of the large ball
        unitBall_samples_int = cell(Mint,1);
        scale_int = zeros(Mint,1);
        for i = 1 : Mint
            scale_int(i) = rand;
            old_dist = GrDist(MU,uniform_samples_int{i},k,d_type)^2;
            unitBall_samples_int{i} = GrGeodesic(MU,...
                uniform_samples_int{i},k,old_dist,...
                big_rad*scale_int(i),d_type); % Defined below.
        end
        data.clean_samples_bb_int = unitBall_samples_int;
        
        %% Generate points from a second uniform distribution
        nSBint = ceil(nSB*int_pts_ratio);
        S2int = randn(n,small_k,nSBint);
        T2int = zeros(size(S2int));
        Q2int = zeros(size(S2int));
        R2int = zeros(small_k,small_k);
        for i = 1 : nSBint
            [Q2int(:,:,i), R2int(:,:,i)] = qr(S2int(:,:,i),0);
            T2int(:,:,i) = Q2int(:,:,i)*diag(sign(diag(R2int(:,:,i))));
        end
        X2int = T2int(:,1:small_k,:);


        % Project these points into the interior of the small ball
        smallBall_samples_int = cell(nSBint,1);
        small_samples_int = cell(nSBint,1);
        scale2_int = zeros(nSBint,1);
        for i = 1:nSBint
            small_samples_int{i} = X2int(:,:,i);
            scale2_int(i) = rand;%small_rad*rand;
            old_dist = GrDist(temp_sc,small_samples_int{i},...
                small_k,d_type)^2;
            smallBall_samples_int{i} = GrGeodesic(temp_sc,...
                small_samples_int{i},small_k,old_dist,...
                small_rad*scale2_int(i),d_type); % Defined below.
        end
        data.clean_samples_sb_int = smallBall_samples_int;
        data.clean_samples_all_int = [data.clean_samples_bb_int;...
            data.clean_samples_sb_int];
        data.clean_samples_both = [data.clean_samples_bb;...
            data.clean_samples_bb_int; data.clean_samples_sb;...
            data.clean_samples_sb_int];
        %% Add a random number of dimensions to each point
        nBBint = nBB*int_pts_ratio; 
        nSBint = nSB*int_pts_ratio; 
        if eqVariance
            S3int = randn(n,kmax,nBBint);
            T3int = zeros(size(S3int));
            Q3int = zeros(size(S3int));
            R3int = zeros(kmax,kmax,nBBint);
            sloppyBall_samples_int = cell(nBBint,1);
            for i = 1 : nBBint
                [oop,~] = qr([unitBall_samples_int{i}, data.center],0);
                [Q3int(:,:,i), R3int(:,:,i)] = qr(S3int(:,:,i)-...
                    oop*(oop'*S3int(:,:,i)),0);
                T3int(:,:,i) = Q3int(:,:,i)*diag(sign(diag(R3int(:,:,i))));
                temp = [unitBall_samples_int{i},...
                    T3int(:,1:dims_bb_int(i)-k,i)] + sqrt(sigmaN)*...
                    (1/(sqrt(dims_bb_int(i))*sqrt(n)))*...
                    randn(n,dims_bb_int(i));
                [temp,~,~] = svd(temp,'econ');
                sloppyBall_samples_int{i} = temp(:,1:dims_bb_int(i));    
            end
            clear dims S3int T3int Q3int R3int
           
            S4int = randn(n,kmax,nSBint);
            T4int = zeros(size(S4int));
            Q4int = zeros(size(S4int));
            R4int = zeros(kmax,kmax,nSBint);
            sloppySmallBall_samples_int = cell(nSBint,1);
            for i = 1 : nSBint
                [oop,~] = qr([smallBall_samples_int{i},...
                    data.small_center, data.center],0);
                [Q4int(:,:,i), R4int(:,:,i)] = qr(S4int(:,:,i)- ...
                    oop*(oop'*S4int(:,:,i)),0);
                T4int(:,:,i) = Q4int(:,:,i)*diag(sign(diag(R4int(:,:,i))));
                temp = [smallBall_samples_int{i}, ...
                    T4int(:,1:max(0,dims_sb_int(i)-small_k),i)] +...
                    sqrt(sigmaN)*((sqrt(max(dims_sb_int(i),small_k))*...
                    1/sqrt(n)))*randn(n,max(dims_sb_int(i),small_k));
                [temp,~,~] = svd(temp,'econ');
                sloppySmallBall_samples_int{i} = temp(:,1:dims_sb_int(i)); 
            end
        else
            T3int = randn(n,kmax,nBBint);
            sloppyBall_samples_int = cell(nBBint,1);
            for i = 1 : nBBint
                T3int(:,:,i) = sqrt(sigmaN)*(1/(sqrt(n)*...
                    sqrt(dims_bb_int(i))))*T3int(:,:,i);
                T3int(:,1:k,i) =  T3int(:,1:k,i) + unitBall_samples_int{i};
                [temp,~,~] = svd(T3(:,1:dims_bb_int(i),i),'econ');
                sloppyBall_samples_int{i} = temp(:,1:dims_bb_int(i));    
            end
            clear dims S3int T3int Q3int R3int

            T4int = randn(n,kmax,nSBint);
            sloppySmallBall_samples_int = cell(nSBint,1);
            for i = 1 : nSBint
                T4int(:,:,i) = sqrt(sigmaN)*(1/(sqrt(n)*...
                    sqrt(dims_sb_int(i))))*T4int(:,:,i);
                T4int(:,1:small_k,i) = T4int(:,1:small_k,i) + ...
                    smallBall_samples_int{i};
                [temp,~,~] = svd(T4(:,1:dims_sb_int(i),i),'econ');
                sloppySmallBall_samples_int{i} = temp(:,1:dims_sb_int(i)); 
            end
        end
        data.noisy_samples_bb_int = sloppyBall_samples_int;
        data.noisy_samples_sb_int = sloppySmallBall_samples_int;
        data.noisy_samples_all_int = [data.noisy_samples_bb_int;...
            data.noisy_samples_sb_int];
        data.noisy_samples_both = [data.noisy_samples_bb; ...
            data.noisy_samples_bb_int; data.noisy_samples_sb; ...
            data.noisy_samples_sb_int];
    end  
end
    


%% HELPER FUNCTIONS:
% Written by T. Marrinan.  Same license applies.

function [out] = SetParam(in)
% ------------------------------------------------------------------------
% SYNTAX:
% [out] = SetParam(in)
%
% DESCRIPTION:
% Checks structure 'in' for required parameter fields and sets the missing
% fields to their default values. Required field is 'default'.
% ------------------------------------------------------------------------
    % Default parameters
    param = struct('default',true,...
        'ambient_dims',100, ...
        'num_samples_bb', 10,...
        'num_samples_sb', 10,...
        'common_dims_bb', 3,...
        'common_dims_sb', 3,...
        'max_dims', 10,...
        'big_rad', .5,...
        'center_sep', (.75)*(.5),...
        'small_rad', (.25)*(.75)*(.5),...
        'noise_var', 0.001,...
        'use_subsample', false,...
        'subsample_percent', .5,...
        'eq_variance', true,...
        'nonuni',false,...
        'nonuni_amt',5,...
        'd_type', 'DPF',...
        'int_pts', false,...
        'int_pts_ratio',1);
%         'use_karcher', false,...
%         'compute_gmeb', false,...

    % read the acceptable names
    paramNames = [fieldnames(param); 'dims_bb'; 'dims_sb'; ...
        'dims_bb_int'; 'dims_sb_int'];
    inNames = fieldnames(in);
    
    % replace defaults with user specified parameters 
    for i = 1 : size(inNames,1)
        inpName = lower(inNames{i});
        if any(strcmp(inpName,paramNames))
            param.(inpName) = in.(inNames{i});
            param.default = false;
        else
            error('%s is not a recognized parameter name',inNames{i})
        end
    end
    
    % check subspace dimensions for compatibility   
    if ~isfield(param,'dims_bb')
        param.dims_bb = randi([param.common_dims_bb,param.max_dims],...
            [param.num_samples_bb,1]); % big ball subspace dimensions
    elseif min(param.dims_bb)<param.common_dims_bb
        warning('min(param.dims_bb) must be larger than param.common_dims_bb. Using random dimensions.');%For random dimensions, remove field with param = rmfield(param,''dims_bb'');')
        param.dims_bb = randi([param.common_dims_bb,param.max_dims],...
            [param.num_samples_bb,1]); % big ball subspace dimensions
    elseif size(param.dims_bb,1) ~= param.num_samples_bb
        warning('param.dims_bb must have size param.num_samples_bb x 1. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_bb'');')
        param.dims_bb = randi([param.common_dims_bb,param.max_dims],...
            [param.num_samples_bb,1]); % big ball subspace dimensions
    elseif max(param.dims_bb>param.max_dims)
        warning('max(param.dims_bb) must be less than param.max_dims. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_bb'');')
        param.dims_bb = randi([param.common_dims_bb,param.max_dims],...
            [param.num_samples_bb,1]); % big ball subspace dimensions
    end
    
    if ~isfield(param,'dims_sb')
        param.dims_sb = randi([param.common_dims_sb,param.max_dims],...
            [param.num_samples_sb,1]); % small ball subspace dimensions
    elseif min(param.dims_sb)<param.common_dims_sb
        warning('min(param.dims_sb) must be larger than param.common_dims_sb. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
        param.dims_sb = randi([param.common_dims_sb,param.max_dims],...
            [param.num_samples_sb,1]); % small ball subspace dimensions
    elseif max(param.dims_sb>param.max_dims)
        warning('max(param.dims_sb) must be less than param.max_dims. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
        param.dims_sb = randi([param.common_dims_sb,param.max_dims],...
            [param.num_samples_sb,1]); % small ball subspace dimensions
    elseif size(param.dims_sb,1) ~= param.num_samples_sb
        warning('param.dims_sb must have size param.num_samples_sb x 1. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
        param.dims_sb = randi([param.common_dims_sb,param.max_dims],...
            [param.num_samples_sb,1]); % small ball subspace dimensions
    end
    if param.int_pts       
        if ~isfield(param,'dims_bb_int')
            if param.nonuni
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil((param.num_samples_bb+param.nonuni_amt)*param.int_pts_ratio),1]); % big ball subspace dimensions
            else
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil(param.num_samples_bb*param.int_pts_ratio),1]); % big ball subspace dimensions
            end
        elseif min(param.dims_bb_int)<param.common_dims_bb
            warning('min(param.dims_bb_int) must be larger than param.common_dims_bb. Using random dimensions.');%For random dimensions, remove field with param = rmfield(param,''dims_bb'');')
            if param.nonuni
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil((param.num_samples_bb+param.nonuni_amt)*param.int_pts_ratio),1]); % big ball subspace dimensions
            else
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil(param.num_samples_bb*param.int_pts_ratio),1]); % big ball subspace dimensions
            end
        elseif size(param.dims_bb_int,1) ~= ceil(param.num_samples_bb*param.int_pts_ratio)
            warning('param.dims_bb_int must have size ceil(param.num_samples_bb*param.int_pts_ratio) x 1. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_bb'');')
            if param.nonuni
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil((param.num_samples_bb+param.nonuni_amt)*param.int_pts_ratio),1]); % big ball subspace dimensions
            else
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil(param.num_samples_bb*param.int_pts_ratio),1]); % big ball subspace dimensions
            end
        elseif max(param.dims_bb_int>param.max_dims)
            warning('max(param.dims_bb_int) must be less than param.max_dims. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_bb'');')
            if param.nonuni
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil((param.num_samples_bb+param.nonuni_amt)*param.int_pts_ratio),1]); % big ball subspace dimensions
            else
                param.dims_bb_int = randi([param.common_dims_bb,param.max_dims],...
                    [ceil(param.num_samples_bb*param.int_pts_ratio),1]); % big ball subspace dimensions
            end
        end

        if ~isfield(param,'dims_sb_int')
            param.dims_sb_int = randi([param.common_dims_sb,param.max_dims],...
                [ceil(param.num_samples_sb*param.int_pts_ratio),1]); % small ball subspace dimensions
        elseif min(param.dims_sb_int)<param.common_dims_sb
            warning('min(param.dims_sb_int) must be larger than param.common_dims_sb. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
            param.dims_sb_int = randi([param.common_dims_sb,param.max_dims],...
                [ceil(param.num_samples_sb*param.int_pts_ratio),1]); % small ball subspace dimensions
        elseif max(param.dims_sb_int>param.max_dims)
            warning('max(param.dims_sb_int) must be less than param.max_dims. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
            param.dims_sb_int = randi([param.common_dims_sb,param.max_dims],...
                [ceil(param.num_samples_sb*param.int_pts_ratio),1]); % small ball subspace dimensions
        elseif size(param.dims_sb_int,1) ~= ceil(param.num_samples_sb*param.int_pts_ratio)
            warning('param.dims_sb_int must have size ceil(param.num_samples_sb*param.int_pts_ratio) x 1. Using random dimensions.');% For random dimensions, remove field param = rmfield(param,''dims_sb'');')
            param.dims_sb_int = randi([param.common_dims_sb,param.max_dims],...
                [ceil(param.num_samples_sb*param.int_pts_ratio),1]); % small ball subspace dimensions
        end
    end
    out = param;
end

function [new_pt] = GrGeodesic(starting_pt, end_pt, k, old_dist, new_dist, metric)
    % ------------------------------------------------------------------------
    % SYNTAX:
    % [new_pt] = GrGeodesic(starting_pt, end_pt, k, old_dist, new_dist)
    %
    % DESCRIPTION:
    % Returns the point 'new_pt' on Gr(k,n) that is 'new_dist' away from 
    % 'starting_point' along the geodesic to 'end_pt'.
    %
    % ------------------------------------------------------------------------
    % REFERENCES:
    % [1] Gallivan, Kyle A., Anuj Srivastava, Xiuwen Liu, and Paul Van Dooren. 
    % "Efficient algorithms for inferences on grassmann manifolds." In IEEE 
    % Workshop on Statistical Signal Processing, 2003, pp. 315-318. IEEE, 2003.
    %
    % NOTES:    
    %
    % LAST EDITED: 22/05/2019 by Tim Marrinan
    % ------------------------------------------------------------------------
    tol = 10^(-6);
    maxiter = 20;
    if strcmp(metric,'DGEO')
        current_dist = GrDist(starting_pt,end_pt,k,metric);
    else
        current_dist = GrDist(starting_pt,end_pt,k,metric)^2;
    end
    alpha_big = 1;
    alpha_little = 0;
    alpha_mid = new_dist/old_dist;
    new_pt = end_pt;
    iter = 0;
    while abs(current_dist - new_dist) > tol && iter < maxiter
        iter = iter + 1;
        %[uc, dc, vc] = svds(starting_pt'*end_pt,k);
        [uc, dc, vc] = svd(starting_pt'*end_pt,'econ');
        uc = uc(:,1:k);
        vc = vc(:,1:k);
        s0 = starting_pt*uc; 
        s1 = end_pt*vc;
        %Theta = real(acos(diag(dc)));
        S = diag(dc);
        T = zeros(k,1); 
        nAngles = min([k,size(S,1)]);
        firstNonzero = k-nAngles;
        for j=1:nAngles
         if S(j)>1-10^(-8)
             if S(j) > 1
                 S(j) = 1;
             end
             T(firstNonzero+j)=sqrt(2*(1-S(j)));    
         else
             T(firstNonzero+j)=real(acos(S(j)));
             %T(k)=acos(S(k));
         end
        end
        Theta = T;

        % This is the suggestion from [1] on how to deal with principal angles
        % equal to zero. It is consistent with the limit of sin(tx)/sin(x) as x ->
        % 0
        omega = zeros(size(Theta));
        z_est = 0.00001;
        omega(Theta>z_est) = sin(alpha_mid*Theta(Theta>z_est))./...
            sin(Theta(Theta>z_est));
        omega(Theta<=z_est) = alpha_mid;

        new_pt = s0*diag(cos(alpha_mid*Theta)) + (s1-s0*diag(cos(Theta)))...
            *diag(omega);
            %*diag(t_inv)*diag(sin(alpha*Theta));
            %*diag(max(0,1./sin(Theta)))*diag(sin(alpha*Theta));
        
        if strcmp(metric,'DGEO')
            current_dist = GrDist(starting_pt,new_pt,k,metric);
        else
            current_dist = GrDist(starting_pt,new_pt,k,metric)^2;
        end
        if current_dist - new_dist > 0
            alpha_big = alpha_mid;
            alpha_mid = (alpha_big+alpha_little)/2;
        else
            alpha_little = alpha_mid;
            alpha_mid = (alpha_big+alpha_little)/2;
        end
    end
end


function [distance]=GrDist(P, Q, N, dist_type)
    % SYNTAX:   [distance]=GrDist(P, Q, N, dist_type)
    %
    % INPUTS:   'P' and 'Q' are matrices that span the two subspaces you wish 
    %           to measure the distance between.  They must be of the same
    %           dimension.
    %
    %           'N' is the number of angles to be used in our distance
    %           calculation. Default value is N = size(P,2);
    %
    %           'dist_type' is a string that specifies which distance measure
    %           you would like to use.  Options are:
    %           [DGEO] for the Geodesic distance
    %           [DPF] for the Projection F-norm
    %           [DC] for the Chordal distance
    %           [DCF] for the Chordal Frobenius
    %           [DFS] for the Fubini-Studi
    %           [MAX] for the maximum angle
    %           [MIN] for the minimum angle
    %           [ALL] for all the angles
    %           [EUC] for the Frobenius norm of P-Q
    %           Default value is 'DGEO'
    %
    % OUTPUTS:  'distance' is the measure between the two subspaces using the
    %           distance specified.
    %
    % USE:      The point of this function is compute the distance between two
    %           subspaces that are thought of as points on a Grassmannian
    %           manifold.
    %
    % NOTES:    I also made it so that I just use as many angles as I can each
    %           time, not a set number.
    %
    %
    %
    % LAST EDITED: 03/13/13 by Tim Marrinan

    % Set up the default values
    distance = 500;

    if ~exist('dist_type','var')
        dist_type = 'DGEO';
    end

    q = min(size(Q,2),size(P,2));
    offset = abs(size(Q,2)-size(P,2));
    N=min(N,q);
    
    if size(P,1)~=size(Q,1)
        disp('GrDist halting due to dimension mismatch.')
        size(P,1)
        size(Q,1)
        return
    end

    %Euclidean
    if strcmp(dist_type,'EUC')
        distance = norm(P-Q,'fro');
        return
    end
    
    % % Test to make sure the matrices are orthonormal
    dd=norm(P'*P-eye(size(P,2)),'fro');
    ee=norm(Q'*Q-eye(size(Q,2)),'fro');

    if dd>1000*eps
        %warning('Basis is not orthonormal. Finding nearest pt. on manifold.')
        %dd
        %msg = 'Error occurred.';
        %error(msg)
        r = size(P,2);
        [P, ~]=qr(P,0);
        P = P(:,1:r);
    end

    if ee>1000*eps
        %warning('Basis is not orthonormal. Finding nearest pt. on manifold.')
        %ee
        %msg = 'Error occurred.';
        %error(msg)
        r = size(Q,2);
        [Q, ~]=qr(Q,0);
        Q = Q(:,1:r);
    end

    % Compute the principle angles between my input matrices and calculate 
    % the distance.
    S = svd(P'*Q);
    T=zeros(N,1);
    %fprintf('\n\tnumber of angles used is %g\n',N);
    for k=1:N
     if S(k)>1-10^(-8)
         if S(k) > 1
             S(k) = 1;
         end
         T(k)=sqrt(2*(1-S(k)));    
     else
         T(k)=real(acos(S(k)));
         %T(k)=acos(S(k));
     end
    end

    %Geodesic distance
    if strcmp(dist_type,'DGEO')
        distance = norm(T);
    end

    %Projection f-norm
    if strcmp(dist_type,'DPF')
        distance = norm(sin(T));
    end

    %chordal distance
    if strcmp(dist_type,'DC')
        distance = norm(sin(T));
    end

    %chordal Frobenius
    if strcmp(dist_type,'DCF')
        distance = norm(2*sin(T/2));
    end

    %Fubini-Studi
    if strcmp(dist_type,'DFS')
        distance = acos(prod(S(1:N)));
    end

    %Min Angle
    if strcmp(dist_type,'MIN')
        distance = T(1);
    end

    %Max Angle
    if strcmp(dist_type,'MAX')
        distance = T(N);
    end

    %All the angles
    if strcmp(dist_type,'ALL')
        distance = T;
    end

    %Max Metric
    if strcmp(dist_type,'MaxMetric')
        if offset>0
            distance = [T;ones(offset,1)*(pi/2)];
        else
            distance = T;
        end
    end
end

