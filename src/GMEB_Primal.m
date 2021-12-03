%   File: GMEB_Primal.m
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
% [iterates, final, options] = GMEB_Primal(data, k, options)
%
% ------------------------------------------------------------------------
% OVERVIEW:
% An implementation of the primal Grassmannian minimum enclosing ball 
% method from [2] with the nonsummable, but square-summable step-size.  
% Cost function is not guarranteed to decrease after each step.
%
% I have no idea how to write the license for my implementation of
% someone else's algorithm.
%
% We wish to solve the primal problem,
%
%       min        max   d(U,X_i)^2
%    U in Gr(p,k)   i 
%
% for some data, {X_i} in Gr(p_i,k) for i=1,...,m.  The distance in this
% problem is the projection Frobenius norm. At optimality, in the support
% set of the center, that is, there should be at least two subspaces that
% achieve the maximum distance to U.  Therefore this primal algorithm
% iterates by taking a geodesic step in the direction of the point furthest
% from the center, with a nonsummable, but square-summable step size. The 
% algorithm terminates after a fixed number of iterations or when two 
% points achieve the the maximum distance.
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
%   'tol'               scalar - threshold for difference in the maximum
%                       distance from the center for the two furthest
%                       poitns
%                       [10^(-4)]
%   'initial_soln'      n x k double - initial value of the primal variable
%                       should be a matrix with orthonormal columns
%                       [ svds([data{:}],k) ] 
%   'maxiter'           scalar - maximum number of outer iterations of 
%                       the algorithm
%                       [500]
%   'per'               scalar - step-size parameter, the numerator in the
%                       step-size determines the length of the first step.
%                       [2]
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
%   distances           nSamp x maxiter double - column t of this matrix
%                       contains the distances between each subspace and
%                       the primal variable at iteration t
%   primal_cost         1 x maxiter double - column t of this matrix
%                       contains the primal cost at iteration t
%   primal_star         1 x maxiter double - column t of this matrix
%                       contains the lowest primal cost at iteration t
%   time                1 x maxiter double - column t of this matrix
%                       contains computation time of iteration t
%   err                 1 x maxiter double - column t of this matrix
%                       contains difference between the largest and second
%                       largest elements of column t of iterates.distances
%
% final     MATLAB structure - contains the best value the variables in the 
%           algorithm when the method terminates
%   Fields:
%   k                   scalar - the subspace dimension of the Grassmannian
%   convergence         boolean - flag indicating whether or not the
%                       duality gap reached the convergence tolerance
%   iterations          scalar - number of iterations the algorithm actual
%                       used before it returned
%   primal_cost         scalar - min(iterates.primal_cost)
%   u                   n x k double - iterates.u_star{final.iterations}
%   time                scalar - sum(iterates.time)
%   err                 scalar - iterates.err(final.iterations)
%
% options   MATLAB structure - contains the parameters for algorithm. 
%           Only differs from input structure if fields were left blank.
%
% -------------------------------------------------------------------------
% DEPENDENCIES:
% [a] SimplexProj.m
%       Not included within. 
%       Written by Nicolas Gillis, see file for documentation.
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
% [3]   Gallivan, Kyle A., Anuj Srivastava, Xiuwen Liu, and Paul Van 
%       Dooren. "Efficient algorithms for inferences on grassmann 
%       manifolds." In IEEE Workshop on Statistical Signal Processing, 
%       2003, pp. 315-318. IEEE, 2003.
%
% ------------------------------------------------------------------------
% CREATED:      13/09/2018 by Tim Marrinan
%
% LAST EDITED:  05/03/2020 by Tim Marrinan
%
% NOTES: 
%
% 04/03/2020 T. Marrinan
% I have no idea how to write the license for my own implementation of
% someone else's algorithm.
%
% ------------------------------------------------------------------------
function [iterates, final, options] = GMEB_Primal(data, k, options)

%% Optional inputs
[nSamp,~] = size(data);
if nargin <= 2
    options = [];
end
if ~isfield(options,'initial_soln') % Initial primal variable
    % Default is equal weights
    [initial_soln,~,~] = svd([data{:}],'econ');
    initial_soln = initial_soln(:,1:k);
    options.initial_soln = initial_soln;
else
    initial_soln = options.initial_soln;
end
if ~isfield(options,'tol') % "Equal" maximum distance threshold
    options.tol = 10^(-4);
    tol = options.tol;
else
    tol = options.tol;
end
if ~isfield(options,'maxiter') % Max iterations of primal algorithm
    options.maxiter = 500;
    maxiter = options.maxiter;
else
    maxiter = options.maxiter;
end
if ~isfield(options,'per') % Step-size parameter
    options.per = 5;
    per = options.per;
else
    per = options.per;
end

%% Check that data matrices are orthonormal (commented for speed tests)
for i = 1 : nSamp
    dd=norm(data{i}'*data{i}-eye(size(data{i},2)));
    if dd>100*eps
        disp('Finding nearest pt. on manifold.')
        [data{i}, ~]=qr(data{i},0);
    end
end

% Initialize output variables
iterates.u = cell(maxiter,1);
iterates.u_star = cell(maxiter,1);
iterates.distances = zeros(nSamp,maxiter);
iterates.primal_cost = zeros(1,maxiter);
iterates.primal_star = zeros(1,maxiter);
iterates.time = zeros(1,maxiter);
iterates.err = zeros(1,maxiter);

final.k = k;
final.convergence = false;
final.iterations = 0;
final.u = 'n/a';
final.primal_cost = 'n/a';
final.time = 0;
final.err = 'n/a';

%% Compute initial solution
t = 1;
tstart = tic;
iterates.u{t} = initial_soln;
% Compute distances
for i = 1 : nSamp
    % Point-to-set projection Frobenius norm, squared
    nEigs = min([k, size(data{i},2)]);
    sig = svd(iterates.u{t}'*data{i},'econ');
    sig = sig(1:nEigs);
    iterates.distances(i,t) = nEigs-sum(sig.^2); %(1/k)*()
end
[dist_sort, ind_sort] = sort(iterates.distances(:,t),'descend');
iterates.u_star{t} = iterates.u{t};
iterates.primal_cost(t) = dist_sort(1);%(ind_sort(1));
iterates.primal_star(t) = iterates.primal_cost(t);
iterates.time(t) = toc(tstart);
iterates.err(t) = dist_sort(1) - dist_sort(2);% tol+1;
imax = ind_sort(1);

%% Main loop
while iterates.err(t) >= tol && t < maxiter
    t = t + 1;
    tstart = tic;
    
    % Take a geodesic step of length delta, according to the method of [2].
    [Uc,Sc,Vc] = svds(iterates.u{t-1}'*data{imax},k); % svds is pretty inefficient sometimes
    S0 = iterates.u{t-1}*Uc;
    S1 = data{imax}*Vc;   
    T = acos(diag(Sc));
    delta = 1/(t+1); % This is the step size
    Gam = diag(cos(T));
    Gam_delta = diag(cos(delta.*T));
    sinverse = zeros(size(T));
    for i = 1:length(T)
        if T(i)~=0
            sinverse(i) = 1./sin(T(i));            
        end
    end
    Sig_inv = diag(sinverse);
    Sig_delta = diag(sin(delta.*T));
    iterates.u{t} = S0*(Gam_delta)+(S1-S0*(Gam))*Sig_inv*Sig_delta;
    
    % Recompute distances    
    for i = 1 : nSamp
        % Point-to-set projection Frobenius norm, squared
        nEigs = min([k, size(data{i},2)]);
        sig = svd(iterates.u{t}'*data{i},'econ');
        sig = sig(1:nEigs);
        iterates.distances(i,t) = nEigs-sum(sig.^2); %(1/k)*()
    end
    [dist_sort, ind_sort] = sort(iterates.distances(:,t),'descend');
    iterates.primal_cost(t) = dist_sort(1);
    [val,b] = min(iterates.primal_cost(1:t));
    iterates.u_star{t} = iterates.u{b};
    iterates.primal_star(t) = val;
    iterates.time(t) = toc(tstart);
    iterates.err(t) = dist_sort(1) - dist_sort(2);% tol+1;
    imax = ind_sort(1);
end

final.iterations = t;
if t < maxiter
    final.convergence = true;
end
final.u = iterates.u_star{t};
final.primal_cost = iterates.primal_star(t);
final.time = sum(iterates.time);
    