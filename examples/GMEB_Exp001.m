%   File: GMEB_Exp001.m
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
% GMEB_Exp001( scenario, nRuns );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This experiment compares the accuracy of the GMEB center computed by the
% proposed algorithm and the existing method of [2] as a function of number
% of iterations and of computation time.
%
% Displays performance plots corresponding to Fig. 5 and Fig. 6 in the
% paper.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig. 5'/'Fig. 6'
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
%       Not included.
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
% [d] GMEB_Primal.m       
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
% [2]   Renard, Emilie, Kyle A. Gallivan, and P-A. Absil. 
%       "A Grassmannian Minimum Enclosing Ball Approach for Common Subspace 
%       Extraction." In International Conference on Latent Variable 
%       Analysis and Signal Separation, pp. 69-78. Springer, Cham, 2018.
%
% ------------------------------------------------------------------------
% CREATED:      03/10/2019 by Tim Marrinan
%
% LAST EDITED:  09/03/2020 by Tim Marrinan
%
% NOTES: 
% ------------------------------------------------------------------------
%clear all; clc;
function GMEB_Exp001( scenario, nRuns )

%% Experimental setup:
%scenario    = 'Fig. 5';                             % Choose scenario
%scenario    = 'Fig. 6';                             % Choose scenario
%nRuns       = 5;                                    % Number of trials
[param, dual_options, primal_options, ~] =  ...
    GMEB_ScenarioSpecification( scenario );         % Set parameters
k           = param.common_dims_bb;                 % Optimal dimension
printout    = 'verbose';                            % 'verbose'/'none'




% Initialize variables
proposed_error  = zeros(nRuns,dual_options.maxiter);
proposed_time   = zeros(nRuns,dual_options.maxiter);
renard_error    = zeros(nRuns,primal_options.maxiter);
renard_time     = zeros(nRuns,primal_options.maxiter);

%% Main loop
for iter = 1 : nRuns
    switch lower(printout)
        case 'verbose'
            clc
            fprintf('\t--------------------------------------------------------\n');
            fprintf('\tExperiment:\t001\n');
            fprintf('\tScenario:\t%s\n',scenario);
            fprintf('\tTrial:\t\t%d/%d\n',iter,nRuns);
            fprintf('\t--------------------------------------------------------\n');
            
            %fprintf('\n\tScenario:\t%s\n',scenario);
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
    [param, dual_options, primal_options, ~] =  ...
        GMEB_ScenarioSpecification( scenario );         % Set parameters

    % Generate data:
    [ data, param ] = GMEB_DataGen( param ) ;

    % Noiseless samples
    pts = data.clean_samples_all;

    % Proposed method
    [proposed_iterates, proposed_optimal, ~] = GMEB_DualSubgrad(pts, k, dual_options);

    % Primal method [1]
    [renard_iterates, renard_optimal, ~] = GMEB_Primal(pts, k, primal_options);

    % Measure the distance from the ground truth center
    for i = 1 : proposed_optimal.iterations
        sig = svd(data.center'*proposed_iterates.u_star{i},0);
        proposed_error(iter,i) = k - sum(sig(1:k).^2);
        proposed_time(iter,i) = sum(proposed_iterates.time(1:i));
    end
    
    % If the algorithm converged without reaching the maximum number of
    % iterations, fix the error for the remaining iterations.
    proposed_error(iter,proposed_optimal.iterations:dual_options.maxiter) = ...
        proposed_error(iter,proposed_optimal.iterations);
    proposed_time(iter,proposed_optimal.iterations:dual_options.maxiter) = ...
        proposed_time(iter,proposed_optimal.iterations);
    
    % Measure the distance from the ground truth center
    for i = 1 : renard_optimal.iterations
        sig = svd(data.center'*renard_iterates.u_star{i},0);
        renard_error(iter,i) = k - sum(sig(1:k).^2);
        renard_time(iter,i) = sum(renard_iterates.time(1:i));
    end
    
    % If the algorithm converged without reaching the maximum number of
    % iterations, fix the error for the remaining iterations.
    renard_error(iter,renard_optimal.iterations:primal_options.maxiter) = ...
        renard_error(iter,renard_optimal.iterations);
    renard_time(iter,renard_optimal.iterations:primal_options.maxiter) = ...
        renard_time(iter,renard_optimal.iterations);           
end


%% Performance plots
MethodColors    = cell(6,1);
MethodColors{1} = [152,78,163]/255;
MethodColors{2} = [141,211,199]/255;
LineStyle{1} = '--';
LineStyle{2} = '-';
l_width = 2;
fig_width = 12;
fig_height = 12;
clear nam
nam{1} = 'Proposed';
nam{2} = '[25]';
error{1} = proposed_error;
error{2} = renard_error;

time{1} = proposed_time;
time{2} = renard_time;
maxiter = min([dual_options.maxiter,primal_options.maxiter]);

% Error vs. Iteration
f{1}    = figure('units','centimeters','Position',...
            [3,3,fig_width,fig_height]);
g = cell(size(nam,2),3);
set(0,'CurrentFigure',f{1})
for i =  size(nam,2):-1:1
    x = 1 : maxiter;
    x2 = [x, fliplr(x)];
    inBetween = [max(error{i}(:,1:maxiter)), fliplr(min(error{i}(:,1:maxiter)))];
    g{i,2} = semilogy(x,median(error{i}(:,1:maxiter)),LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);
    hold on
    g{i,1} = fill(x2, inBetween, MethodColors{i},'edgecolor','none');
    hold on
    set(g{i,1},'facealpha',.3)

end
legend([g{1,2},g{2,2}],nam{1},nam{2});
figtitle = strcat(scenario,'a');
title(figtitle)
ylabel('Error')
xlabel('Iteration')
hold off
               
% Error vs. Cumulative Computation Time
f{2}    = figure('units','centimeters','Position',...
    [3,3,fig_width,fig_height]);

end_time = min([time{1}(:,dual_options.maxiter);time{2}(:,primal_options.maxiter)]);
spots1 = time{1}<end_time;
[~,ind1] = min(spots1,[],2);
spots2 = time{2}<end_time;
[~,ind2] = min(spots2,[],2);
steps = min([ind1;ind2]);
t_int = linspace(0,end_time,steps);
timeError = cell(2,1);
for i = 1 : 2
    timeError{i} = ones(nRuns,steps);
    for j = 1 : nRuns
        timeError{i}(j,:) = max(error{i}(j,:));
        tIter = zeros(1,steps);
        for l = 1 : steps
            if sum(time{i}(j,:)<=t_int(l))>0
                tIter(l) = find(time{i}(j,:)<=t_int(l),1,'last');
                 timeError{i}(j,l) = error{i}(j,tIter(l));
            end
        end

    end
end     

w = cell(size(nam,2),3);
set(0,'CurrentFigure',f{2})
for i =  1:size(nam,2)
    x = t_int;
    x2 = [x, fliplr(x)];
    inBetween = [min(timeError{i}), fliplr(max(timeError{i}))];
    w{i,2} = semilogy(x,median(timeError{i}),LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);
    hold on
    w{i,1} = fill(x2, inBetween, MethodColors{i},'edgecolor','none');
    hold on
    set(w{i,1},'facealpha',.3)
end
legend([w{1,2},w{2,2}],nam{1},nam{2});
figtitle = strcat(scenario,'b');
title(figtitle)
ylabel('Error')
xlabel('Cumulative time')
xlim([0,end_time])
hold off