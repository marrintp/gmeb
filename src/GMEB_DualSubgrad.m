%   File: GMEB_DualSubgrad.m
%   Copyright (c) <2020> <University of Mons>
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
%       published under this license. The University of Mons
%       is granted the non-exclusive right to publish modifications
%       or contributions in future versions of the Software free of charge.
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
%
% ------------------------------------------------------------------------
% SYNTAX:  
% [iterates, final, options] = GMEB_DualSubgrad(data, k, options);
%
% ------------------------------------------------------------------------
% OVERVIEW: 
% We wish to solve the primal problem,
%
%       min        max   d(U,X_i)^2
%    U in Gr(p,k)   i 
%
% for some data, {X_i} in Gr(p_i,k) for i=1,...,m.  The distance in this
% problem is the projection Frobenius norm, and the optimal solution, U*, 
% can be computed as the dominant k-dimensional eigenspace of the weighted 
% sum, 
%
%       sum lambda_i X_i X_i^T  =   U D U^T.
%        i
% This would be easy to compute, except the weights, lambda_i, are unknown.
% This code solves the Lagrange dual of this problem, which can be written
% as
%
%       max         q(lambda)
%   lambda in R^m
%
%       s.t.        q(lambda)    =   min        sum ( lambda_i d(U,X_i)^2 )
%                                U in Gr(p,k)    i
%                   lambda_i    >=   0 for i = 1,...,m
%                || lambda ||_1  =   1.
% 
% This dual can be solved via a 'subgradient' algorithm, where the scare
% quotes are there because we are actually looking at a supergradient since
% this is a maximization problem.  For a fixed lambda, let
%
%   g_i(U) = lambda_i d^2(U,X_i)
%
% a vector included in the 'subdifferential' is then the vector of function
% values at the minimum of the sum.  That is 
%
%   g^{(t)} = [g_1(U*), g_2(U*), ..., g_m(U*)]^T in del q(lambda).
%
% The subgradient algorithm in general terms is then:
% (1) Intialize the weights
%       lambda^(0) = unit vector with equal elements
%       j = 0
%
% Repeat until convergence:
%
% (2) Solve (via the eigenvalue decomposition mentioned above)
%       min        sum ( lambda_i^(j) d^2(U,X_i) )
%   U in Gr(p,k)    i
%
% (3) Pick a subgradient
%        g^(j) = [f_1^(j)(U*), f_2^(j)(U*), ..., f_m^(j)(U*)]^T
%
% (4) Take a step in the direction of the subgradient with stepsize alpha
% (to be decided)
%       \hat{lambda}^(j+1) = lambda^(j) - alpha del g^(j)
%
% (5) Project onto the feasible set (using Nicolas' simplex projection
% function)
%       lambda^(j+1) = P( \hat{lambda^(j+1)} )
%
% ------------------------------------------------------------------------
% INPUTS:
% data      M x 1 cell array of subspace bases - Cells contain matrices of 
%           size n x p_i that have orthonormal columns.
%
% k         scalar - subspace dimension of the Grassmannian on which to
%           compute the minimum enclosing ball, Gr(k,n)
%
% options   MATLAB structure - contains the parameters for algorithm. 
%           If not included, default parameters will be used.
%   Fields and [default values]:
%   'tol'               scalar - convergence threshold for the duality gap
%                       and for lack of movement of the dual variable 
%                       [10^(-4)]
%   'step_type'         string - chooses the step size. options are 
%                       'nonSum-sqSum' for per/t or 'nonSum-diminish' for
%                       per/sqrt(t).
%                       ['nonSum-sqSum']
%   'useLineSearch'     boolean - line-search flag. if true, algorithm uses
%                       a backtracking line-search in an attempt to find an
%                       iterate that increases the dual cost
%                       [true]
%   'initial_soln'      nSamp x 1 double - initial value of the dual 
%                       variable
%                       [(1/nSamp).*(ones(nSamp,1))]
%   'maxiter'           scalar - maximum number of outer iterations of 
%                       subgradient algorithm
%                       [200]
%   'per'               scalar - step-size parameter, the numerator in the
%                       step-size determines the length of the first step.
%                       [2]
%   'growth'            scalar - growth parameter for line-search.
%                       determines the amount that the step-size increases
%                       after a step is accepted in the back-tracking 
%                       line-search. should be in the range (1,2).
%                       [1.9]
% 'findStationaryPoint' boolean - flag that returns solution if dual 
%                       variables don't change in 10 iterations.
%                       [false]
%
% ------------------------------------------------------------------------
% OUTPUTS:
% iterates  MATLAB structure - contains the value the variables in the dual
%           subgradient algorithm at each iteration.
%   Fields:
%   u                   maxiter x 1 cell array - cell t contains the 
%                       primal variable at iteration t, an n x k matrix. 
%   u_star              maxiter x 1 cell array - cell t contains the 
%                       primal variable that contains the lowest primal
%                       cost at iteration t, an n x k matrix.
%   lambda              nSamp x maxiter double - column t of this matrix
%                       contains the dual variable at iteration t
%   lambda_star         nSamp x maxiter double - column t of this matrix
%                       contains the dual variable with the highest dual 
%                       cost at iteration t
%   g                   nSamp x maxiter double - column t of this matrix
%                       contains the subgradient at iteration t
%   alpha               1 x maxiter double - column t of this matrix
%                       contains the accepted step length at iteration t
%   dual_cost           1 x maxiter double - column t of this matrix
%                       contains the dual cost at iteration t
%   primal_cost         1 x maxiter double - column t of this matrix
%                       contains the primal cost at iteration t
%   dual_star           1 x maxiter double - column t of this matrix
%                       contains highest dual cost at iteration t
%   primal_star         1 x maxiter double - column t of this matrix
%                       contains the lowest primal cost at iteration t
%   duality_gap         1 x maxiter double - column t of this matrix
%                       contains dual_star(t) - primal_star(t)
%   time                1 x maxiter double - column t of this matrix
%                       contains computation time of iteration t
%   movement            1 x maxiter double - column t of this matrix
%                       contains dual_star(t) -dual_star(t-10)
%
% final     MATLAB structure - contains the best value the variables in the 
%           dual subgradient algorithm when the algorithm terminates
%   Fields:
%   k                   scalar - the subspace dimension of the Grassmannian
%   convergence         boolean - flag indicating whether or not the
%                       duality gap reached the convergence tolerance
%   iterations          scalar - number of iterations the algorithm actual
%                       used before it returned
%   dual_cost           scalar - max(iterates.dual_cost)
%   primal_cost         scalar - min(iterates.primal_cost)
%   duality_gap         scalar - final.dual_cost - final.primal_cost
%   u                   n x k double - iterates.u_star{final.iterations}
%   lambda              nSamp x 1 double - 
%                       iterates.lambda_star(:,final.iterations)
%   time                scalar - sum(iterates.time)
%
% options   MATLAB structure - contains the parameters for algorithm. 
%           Only differs from input structure if fields were left blank.
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [a] SimplexProj.m
%       Not included within. 
%       Written by Nicolas Gillis. 
%       See file for documentation and reference [2] for details.
%       
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% [2]   Renard, Emilie, Kyle A. Gallivan, and P-A. Absil. 
%       "A Grassmannian Minimum Enclosing Ball Approach for Common 
%       Subspace Extraction." In International Conference on Latent 
%       Variable Analysis and Signal Separation, pp. 69-78. Springer, 
%       Cham, 2018.
%
% [3]   Gillis, Nicolas. 
%       "Successive nonnegative projection algorithm for robust 
%       nonnegative blind source separation." SIAM Journal on Imaging 
%       Sciences 7, no. 2 (2014): 1420-1450.
%
% ------------------------------------------------------------------------
% CREATED:      21/02/2018 by Tim Marrinan
%
% LAST EDITED:  05/03/2020 by Tim Marrinan
%
% NOTES: 
%
% 04/03/2020 T. Marrinan
% Currently the code is the negative of what is written in the paper.  That
% is, it solves a primal minimization and a dual maximization problem.
% This is really what we want, but the paper has it the other way around so
% it can use the term "subgradient" accurately.
%
% Should rewrite code to match paper, but I don't really want to.
%
% ------------------------------------------------------------------------

function [iterates, final, options] = GMEB_DualSubgrad(data, k, options)
% assumes data is a cell array of size nSamp x 1 and each cell
% is an orthonormal matrix of size nSamples x nDims

%% Optional inputs
[nSamp,~] = size(data);
if nargin <= 2
    options = [];
end
if ~isfield(options,'tol') % Duality gap threshold
    options.tol = 10^(-4);
    tol = options.tol;
else
    tol = options.tol;
end
if ~isfield(options,'step_type') % Subgradient step type
    % Default is a nonsummable, but square-summable step
    options.step_type = 'nonSum-sqSum'; % 'nonSum-diminish';%
    step_type = options.step_type;
else
    step_type = options.step_type;
end
if ~isfield(options,'useLineSearch') % Line-search flag
    options.useLineSearch = false;
    useLineSearch = options.useLineSearch;
else
    useLineSearch = options.useLineSearch;
end
if ~isfield(options,'initial_soln') % Initial dual variable
    % Default is equal weights
    options.initial_soln = (1/nSamp).*(ones(nSamp,1));
    initial_soln = options.initial_soln;
else
    initial_soln = options.initial_soln;
end
if ~isfield(options,'maxiter') % Max iterations of subgrad algorithm
    options.maxiter = 200;
    maxiter = options.maxiter;
else
    maxiter = options.maxiter;
end
if ~isfield(options,'per') % Step-size parameter
    options.per = 2;
    per = options.per;
else
    per = options.per;
end
if ~isfield(options,'growth') % Growth parameter for line-search
    % Increases step-size if a step is accepted
    options.growth = 1.9;
    growth = options.growth;
else
    growth = options.growth;
end
if ~isfield(options,'findStationaryPoint') % Return solution if iterates don't change
    options.findStationaryPoint = false;
    findStationaryPoint = options.findStationaryPoint;
else
    findStationaryPoint = options.findStationaryPoint;
end


%% Check that data matrices are orthonormal (commented for speed tests)
for i = 1 : nSamp
    dd=norm(data{i}'*data{i}-eye(size(data{i},2)));
    if dd>100*eps
        disp('Finding nearest pt. on manifold.')
        [data{i}, ~]=qr(data{i},0);
    end
end


%% Initialize variables
iterates.u = cell(maxiter,1);
iterates.u_star = cell(maxiter,1);
iterates.lambda = zeros(nSamp,maxiter);
iterates.lambda_star = zeros(nSamp,maxiter);
iterates.g = zeros(nSamp,maxiter);
iterates.alpha = zeros(1,maxiter);
iterates.dual_cost = zeros(1,maxiter);
iterates.primal_cost = zeros(1,maxiter);
iterates.primal_star = zeros(1,maxiter);
iterates.duality_gap = zeros(1,maxiter);
iterates.dual_star = zeros(1,maxiter);
iterates.time = zeros(1,maxiter);
iterates.movement = zeros(1,maxiter);


final.k = k;
final.convergence = false;
final.duality_gap = 'n/a';
final.primal_cost = 'n/a';
final.dual_cost = 'n/a';
final.iterations = 0;
final.u = 'n/a';
final.lambda = 'n/a';
final.time = 0;


%% Compute initial solution
t = 1;
tstart = tic;
lambda = initial_soln;  % Dual variable
support = lambda>0; % Nonzero dual variables
result = cellfun(@times,data(support),num2cell(sqrt(lambda(support))),'uni',false); % Weighted bases
concatenatedBases =[result{:}];
%[u,~,~] = svd(concatenatedBases,'econ'); % Primal variable
clear u
[u,~,~] = svd(concatenatedBases,'econ'); % primal variable
if k > size(u,2)
    [v,~,~] = svd( [data{:}]-u*(u'*[data{:}]) ,'econ');
    u = [u, v(:,1:k-size(u,2))];
else
    u = u(:,1:k);
end
g = zeros(nSamp,1);
% Compute distances
for i = 1 : nSamp
    % Point-to-set projection Frobenius norm, squared
    nEigs = min([k, size(data{i},2)]); 
    sig = svd(u'*data{i},'econ');
    sig = sig(1:nEigs);
    g(i) = nEigs-sum(sig.^2); % g is the negative of g from the paper %(1/k)*
end
alpha = per/sqrt(t); % Step size

% Store data
iterates.u{t} = u;
iterates.u_star{t} = u;
iterates.lambda(:,t) = lambda;
iterates.lambda_star(:,t) = lambda;
iterates.g(:,t) = g;
iterates.alpha(t) = alpha;
iterates.dual_cost(t) = lambda'*g;
iterates.dual_star(t) = iterates.dual_cost(t);
iterates.primal_cost(t) = max(g);
iterates.primal_star(t) = iterates.primal_cost(t);
iterates.duality_gap(t) = iterates.primal_star(t)-iterates.dual_star(t);
iterates.time(t) = toc(tstart);
iterates.movement(1:10) = 1;

final.u = u;
final.lambda = lambda;
final.primal_cost = max(g);
final.dual_cost = iterates.dual_cost(t);
final.duality_gap = iterates.duality_gap(t);
final.time = iterates.time(t);

%% Main loop: stop when duality gap = 0, step size = 0, or iter = max_iters
while iterates.duality_gap(t) > tol  && t < maxiter 
    tstart = tic;
    t = t + 1;
    switch lower(step_type)
        case 'nonsum-sqsum'
            alpha = per/t;
        case 'nonsum-diminish'
            alpha = per/sqrt(t);
        otherwise
            error('Unknown step_type'); 
    end

    % Take subgradient step
    lambda_temp = SimplexProj(lambda + alpha*g); % Update dual variables
    support = lambda_temp>0; % Support set = nonzero dual variables
    result = cellfun(@times,data(support),num2cell(sqrt(lambda_temp(support))),...
        'uni',false);
    concatenatedBases =[result{:}];
    clear u
    [u,~,~] = svd(concatenatedBases,'econ'); % Update primal variable
    if k > size(u,2)
        [v,~,~] = svd( [data{:}]-u*(u'*[data{:}]) ,'econ');
        u = [u, v(:,1:k-size(u,2))];
    else
        u = u(:,1:k);
    end
    g_temp = zeros(nSamp,1);
    for i = 1 : nSamp % Update distances                    
        nEigs = min([k, size(data{i},2)]);
        sig = svd(u'*data{i},'econ');
        sig = sig(1:nEigs);
        g_temp(i) = nEigs-sum(sig.^2); % (1/k)*()
    end

    % Take a shorter step if the dual cost does not increase 
    % (should probably update this to take a longer step if possible)
    if useLineSearch
        dual_cost_temp = lambda_temp'*g_temp;
        while dual_cost_temp <= iterates.dual_cost(t-1) && alpha > 10^(-4)*iterates.alpha(t)
            %alpha = alpha/2; % Cut the step size in half
            per = per/2; % Cut the step size in half
            switch lower(step_type)
                case 'nonsum-sqsum'
                    alpha = per/t;
                case 'nonsum-diminish'
                    alpha = per/sqrt(t);
                otherwise
                    error('Unknown step_type'); 
            end
            lambda_temp = SimplexProj(lambda + alpha*g); % g is the negative of g in the paper.
            % Check for increased dual cost
            support = lambda_temp>0;
            result = cellfun(@times,data(support),...
                num2cell(sqrt(lambda_temp(support))),'uni',false);
            concatenatedBases =[result{:}];
            clear u
            [u,~,~] = svd(concatenatedBases,'econ'); % Update primal variable
            if k > size(u,2)
                [v,~,~] = svd( [data{:}]-u*(u'*[data{:}]) ,'econ');
                u = [u, v(:,1:k-size(u,2))];
            else
                u = u(:,1:k);
            end
            g_temp = zeros(nSamp,1);
            for i = 1 : nSamp
                nEigs = min([k, size(data{i},2)]);
                sig = svd(u'*data{i},'econ');
                sig = sig(1:nEigs);
                g_temp(i) =nEigs-sum(sig.^2); %(1/k)*()
            end
            dual_cost_temp = lambda_temp'*g_temp;
        end
        % I actually don't know what the point of this is.  
        % Since the step size is independent of anything relating to 
        % the subgradient, this doesn't do anything.
        %alpha = growth*alpha;
        per = growth*per; 
        lambda = lambda_temp;
        g = g_temp;
    else
        lambda = lambda_temp;
        g = g_temp;
    end

    % Store updated variables
    iterates.u{t} = u;
    iterates.lambda(:,t) = lambda;
    iterates.g(:,t) = g;
    iterates.alpha(t) = alpha;
    iterates.dual_cost(t) = lambda'*g;
    iterates.primal_cost(t) = max(g);
    [iterates.dual_star(t),ind] = max(iterates.dual_cost(1:t));
    iterates.lambda_star(:,t) = iterates.lambda(:,ind);
    [val,b] = min(iterates.primal_cost(1:t));
    iterates.primal_star(t) = val;
    iterates.u_star{t} = iterates.u{b};
    iterates.duality_gap(t) =  iterates.primal_star(t) ...
        - iterates.dual_star(t);
    iterates.time(t) = toc(tstart);

    if t > 10
        iterates.movement(t) = abs(iterates.dual_star(t-10) - iterates.dual_star(t));
    end
    if findStationaryPoint
        if iterates.movement(t) <= tol
            final.u = iterates.u_star{t};
            final.lambda = iterates.lambda_star(:,t);
            final.primal_cost = iterates.primal_star(t);
            final.dual_cost = iterates.dual_star(t);
            final.duality_gap = iterates.duality_gap(t);
            if final.duality_gap <= tol
                final.convergence = true;
            end
            final.iterations = t;    
            final.time = sum(iterates.time);
            return;
        end
    end
end


final.u = iterates.u_star{t};
final.lambda = iterates.lambda_star(:,t);
final.primal_cost = iterates.primal_star(t);
final.dual_cost = iterates.dual_star(t);
final.duality_gap = iterates.duality_gap(t);
if final.duality_gap <= tol
    final.convergence = true;
end
final.iterations = t;    
final.time = sum(iterates.time);