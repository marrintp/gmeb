%   File: GMEB_OrderSelection.m
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
% [penalty, hybrid, options] = GMEB_OrderSelection(data, options)
%
% ------------------------------------------------------------------------
% OVERVIEW: 
% This function selects the order of the common subspace in a collection of
% linear subspace data based on the center of the minimum enclosing ball.
% There are two order-selection rules implemented by this function.  
%
% 'Penalty' method:
% Minimizes the sum of the value of the objective function from the 
% minimum enclosing ball problem and the value of a penalty term that 
% measures the amount of common information not included in the central 
% subspace. If U*(k) is the center of the minimum enclosing ball on 
% Gr(k,n), and the d( , ) is the squared chordal distance, the two terms 
% in the sum are:
%
%           0 for k = 0 
% cObj(k) =             
%           max_{i=1,...,M} (1/k) d( U*(k), X_i ) for k=1,...,max{dim(X_i)}
%
%
%
%           1 for k = 0
% cPen(k) =             
%           min_{j=1,...,M} 1 - (1/min{n-k,dim(X_j)}) d(U*_{perp}(k),X_j)
%           for k = 1,..., max{dim(X_j)},
%
% and then the order selection rule is
%
%         k* = argmin_{k=0,...max{dim(X_i)}} cObj(k) + cPen(k).
%
% 'Hybrid' method:
% Combines the eigenvalue threshold of [1] and [2] with the minimum
% enclosing ball problem of [3]. Let lambda*(k) be a vector of dual
% variables such that the dominant k eigenvectors of 
%
%                      Sum lambda_i*(k) X_i X_i^T
%
% form the center of the minimum enclosing ball on Gr(k,n) and let 
%
%           d_1(k) >= d_2(k) >= ... >= d_R(k) be the eigenvalues.
%
% (For k=0 let lambda_i(0) = 1/M for all i). Define
%
%           sum_{r=k+1,...,R} d_r(k) for k = 0
% E(k) = 
%           sum_{r=1,...,k} 1-d_r(k) + sum_{r=k+1,...,R} d_r(k) for k =
%           1,...,max{dim(X_i)}.
%
% Then the hybrid order selection rule is
%
%       k* = argmin_{k=0,...max{dim(X_i)}} E(k).
%
% ------------------------------------------------------------------------
% INPUTS:
% data      M x 1 cell array of subspace bases - Cells contain matrices of 
%           size n x p_i that have orthonormal columns.
%
% options   MATLAB structure - contains the parameters for GMEB_DualSubgrad 
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
% penalty   MATLAB structure - contains the parameters for GMEB_DualSubgrad 
%           If not included, default parameters will be used.
%   Fields:
%   'kStar'             scalar - the selected order of the GMEB center
%                       using the penalty rule (index of the minimum of 
%                       penalty.cost) - 1
%   'cStar'             scalar - the cost of the selected order using the
%                       penalty rule. minimum of penalty.cost
%   'cObj'              kmax+1 x 1 double - row t contains the value of 
%                       the objective term for k = t-1
%   'cPen'              kmax+1 x 1 double - row t contains the value of 
%                       the penalty term for k = t-1
%   'cost'              kmax+1 x 1 double - penalty.cObj + penalty.cPen
%
% hybrid    MATLAB structure - contains the parameters for GMEB_DualSubgrad 
%           If not included, default parameters will be used.
%   Fields:
%   'kStar'             scalar - the selected order of the GMEB center
%                       using the hybrid rule (index of the minimum of 
%                       penalty.cost) - 1
%   'cStar'             scalar - the cost of the selected order using the
%                       hybrid rule. minimum of hybrid.cost
%   'cObj'              kmax+1 x 1 double - row t contains the value of 
%                       the objective term for k = t-1
%   'cPen'              kmax+1 x 1 double - row t contains the value of 
%                       the penalty term for k = t-1
%   'cost'              kmax+1 x 1 double - penalty.cObj + penalty.cPen
%
% options   MATLAB structure - contains the parameters for GMEB_DualSubgrad 
%           Only differs from input structure if fields were left blank.
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [a] GMEB_DualSubgrad.m
%       Not included. 
%       Written by Tim Marrinan, see file for documentation.
%       
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% [2]   Garg, Vaibhav, Ignacio Santamaria, David Ramirez, and Louis Scharf. 
%       "Subspace Averaging and Order Determination for Source 
%       Enumeration." IEEE Transactions on Signal Processing (2019).
%
% [3]   Santamaría, Ignacio, Louis L. Scharf, Chris Peterson, Michael 
%       Kirby, and J. Francos. "An order fitting rule for optimal subspace 
%       averaging."  In 2016 IEEE Statistical Signal Processing Workshop 
%       (SSP), pp. 1-4. IEEE, 2016.
%
% ------------------------------------------------------------------------
% CREATED:      05/03/2020 by Tim Marrinan
%
% LAST EDITED:  05/03/2020 by Tim Marrinan
%
% NOTES: 
%
% 05/03/2020 T. Marrinan 
% This used to be part of the GMEB code, but I made it a stand-alone
% function for timing and modularity.
% ------------------------------------------------------------------------

function [penalty, hybrid, options] = GMEB_OrderSelection(data, options)
%% Optional inputs
% Options can also be specified for 'GMEB_DualSubgrad.m' or the default
% parameters will be used.
[nSamp,~] = size(data);
if nargin <= 1
    options = [];
end
% Maximum dimension to check
dims = zeros(nSamp,1);
for i = 1 : nSamp
    dims(i) = size(data{i},2);
end
kmax = max(dims);

% Initialize variables
penalty.kStar = 'n/a';
penalty.cStar = 'n/a';
penalty.cObj = zeros(kmax+1,1);
penalty.cPen = zeros(kmax+1,1);
penalty.cost = zeros(kmax+1,1);

hybrid.kStar = 'n/a';
hybrid.cStar = 'n/a';
hybrid.cObj = zeros(kmax+1,1);
hybrid.cPen = zeros(kmax+1,1);
hybrid.cost = zeros(kmax+1,1);

[uTot,~,~] = svd([data{:}],'econ');
dimTot = size(uTot,2);
oVec = zeros(nSamp,kmax);
pVec = zeros(nSamp,kmax);        

for k = 1 : kmax
    [~, optimal, ~] = GMEB_DualSubgrad(data, k, options);
    
    % Penalty method
    center = optimal.u;
    [uPerp,~,~] = svd(uTot-center*(center'*uTot),'econ');
    for i = 1 : nSamp
        % Compute objective terms
        nEigs = min([k, size(data{i},2)]);
        sig = svd(center'*data{i},'econ');
        sig = sig(1:nEigs);
        oVec(i,k) =  (1/k)*(nEigs-sum(sig.^2));
        % Compute penalty terms
        if dimTot-k > 0            
            nEigs = min([dimTot-k, size(data{i},2)]);
            sig = svd(uPerp(:,1:dimTot-k)'*data{i},'econ');
            sig = sig(1:nEigs);
            pVec(i,k) = (1/nEigs)*(nEigs-sum(sig.^2));
        end
    end
    if k == 1
        penalty.cObj(k) = 0;
        penalty.cPen(k) = 1;
       penalty.cost(k) = penalty.cObj(k) + penalty.cPen(k); 
    end
    penalty.cObj(k+1) = max(oVec(:,k));
    penalty.cPen(k+1) = 1-max(pVec(:,k));
    penalty.cost(k+1) = penalty.cObj(k+1) + penalty.cPen(k+1);

    % Hybrid method
    support = optimal.lambda>0;
    result = cellfun(@times,data(support),...
        num2cell(sqrt(optimal.lambda(support))),'uni',false);
    sig = svd([result{:}],'econ');
    sig = (sig.^2);
    if k == 1
        hybrid.cObj(k) = 0;
        hybrid.cPen(k) = sum(sig); 
        hybrid.cost(k)= hybrid.cObj(k) + hybrid.cPen(k);       
    end
    hybrid.cObj(k+1) = sum(1-sig(1:k));
    hybrid.cPen(k+1) = sum(sig(k+1:end));
    hybrid.cost(k+1) = hybrid.cObj(k+1) + hybrid.cPen(k+1); 
end
[penalty.cStar,ind] = min(penalty.cost);
penalty.kStar = ind - 1;

[hybrid.cStar,ind] = min(hybrid.cost);
hybrid.kStar = ind - 1;