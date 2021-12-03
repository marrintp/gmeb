%   File: GMEB_Exp002.m
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
% ------------------------------------------------------------------------
% SYNTAX:
% GMEB_Exp002( scenario, nRuns );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This experiment compares the number of iterations needed to find a 
% stationary point of the dual subgradient algorithm with a naive 
% initialization vs. a warm-start.
%
% Displays performance plots corresponding to Figure 7.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig. 7a'/'Fig. 7b'
%
% nRuns     Scalar that specficies the number of Monte Carlo trials to run.
%           (Paper uses 100)
%
% ------------------------------------------------------------------------
% OUTPUTS: 
% n/a
%
% ------------------------------------------------------------------------
% DEPENDENCIES:
% [a] GMEB_ScenarioSpecification.m
%       Not included within.
%       Sets scenario parameters.
%       Written by T. Marrinan; see file for documentation.
%
% [b] GMEB_DataGen.m
%       Not included.
%       Written by T. Marrinan, see file for documentation.
%
% [c] GMEB_DualSubgrad.m
%       Not included. 
%       Written by Tim Marrinan, see file for documentation.
%
% [d] aboxplot.m       
%       Not included.
%       Written by Alex Bikfalvi, see file for documentation.
%       http://alex.bikfalvi.com/research/advanced_matlab_boxplot/
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% ------------------------------------------------------------------------
% CREATED:      21/10/2019 by Tim Marrinan
%
% LAST EDITED:  09/03/2030 by Tim Marrinan
%
% NOTES: 
% ------------------------------------------------------------------------
%clear all; clc;
function GMEB_Exp002( scenario, nRuns )

%% Experimental setup:
%scenario    = 'Fig. 7a';                            % Choose scenario
%scenario    = 'Fig. 7b';                            % Choose scenario
%nRuns       = 100;                                  % Number of trials 
[param, dual_options, ~, ~] =  ...
    GMEB_ScenarioSpecification( scenario );         % Set parameters
kmax        = param.max_dims;                       % Max dimension
printout    = 'verbose';                            % 'verbose'/'none'


% Initialize variables
naive_iterations = zeros(nRuns,param.max_dims);
naive_convergence = zeros(nRuns,param.max_dims);
naive_time = zeros(nRuns,param.max_dims);
warm_iterations = zeros(nRuns,param.max_dims);
warm_convergence = zeros(nRuns,param.max_dims);
warm_time = zeros(nRuns,param.max_dims);

%% Main loop
for iter = 1 : nRuns
    switch lower(printout)
        case 'verbose'
            clc
            fprintf('\t--------------------------------------------------------\n');
            fprintf('\tExperiment:\t002\n');
            fprintf('\tScenario:\t%s\n',scenario);
            fprintf('\tTrial:\t\t%d/%d\n',iter,nRuns);
            fprintf('\t--------------------------------------------------------\n');
        otherwise
    end

    % Need to clear the subspace dimensions from the parameter
    % vector each time otherwise they will not be randomly
    % resampled
    if isfield(param,'dims_bb')
        param = rmfield(param,'dims_bb');
    end
    if isfield(param,'dims_sb')
        param = rmfield(param,'dims_sb');
    end
    if isfield(param,'dims_bb_int')
        param = rmfield(param,'dims_bb_int');
    end
    if isfield(param,'dims_sb_int')
        param = rmfield(param,'dims_sb_int');
    end
    
    % Reset the parameters each time
    % (The nonuniform option is a bit wonky and adds more points each time.
    % I should probably fix that, huh?)
    [param, dual_options, ~, ~] =  ...
        GMEB_ScenarioSpecification( scenario );         % Set parameters
    
    % Generate data:
    [ data, param ] = GMEB_DataGen( param ) ;

    % Noisy samples:
    % In this case the central subspace has no noise, but this allows us to
    % add the random other dimenisons.
    pts = data.noisy_samples_all;

    % Solve the GMEB problem for all values of k
    n_sets = size(pts,1);
    initial_soln = (1/n_sets).*(ones(n_sets,1));
    initial_soln_old = (1/n_sets).*(ones(n_sets,1));
    for i = 1:kmax
        dual_options.initial_soln = initial_soln;
        [~, optimal, ~] = GMEB_DualSubgrad(pts, i, dual_options);
        naive_iterations(iter,i) = optimal.iterations;
        naive_convergence(iter,i) = optimal.convergence;
        naive_time(iter,i) = optimal.time;

        dual_options.initial_soln = initial_soln_old;
        [~, optimal, ~] = GMEB_DualSubgrad(pts, i, dual_options);
        initial_soln_old = optimal.lambda;
        warm_iterations(iter,i) = optimal.iterations;
        warm_convergence(iter,i) = optimal.convergence;
        warm_time(iter,i) = optimal.time;
    end

end


%% Performance plots
MethodColors    = cell(2,1);
MethodColors{1} = [253,180,98]/255;
MethodColors{2} = [228,26,28]/255;

fig_width = 12;
fig_height = 12;

f{1}    = figure('units','centimeters','Position',...
            [3,3,fig_width,fig_height]);

clear nam
nam{1} = 'Naive Initialization';
nam{2} = 'Warm Start';
iterations{1} = naive_iterations(:,2:kmax);
iterations{2} = warm_iterations(:,2:kmax);


set(0,'CurrentFigure',f{1})
aboxplot(iterations,'labels',2:kmax,'colormap',[MethodColors{1};MethodColors{2}])
legend(nam{1},nam{2})
title(scenario)
ylabel('Iterations')
xlabel('k')
ylim([0,60])
