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
% GMEB_Exp003( scenario, nRuns );
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This experiment compares different methods for estimating the optimal 
% dimension for averaging subspace data.
%
% Displays performance plots corresponding to Figure 8-10.
%
% ------------------------------------------------------------------------
% INPUTS:
% scenario  String that specificies the figure to reproduce.
%           'Fig. 8'/'Fig. 9'/'Fig. 10'
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
% [d] knee_pt.m       
%       Not included.
%       Written by Dmitry Kaplan, see file for documentation.
%       Dmitry Kaplan (2020). Knee Point (https://www.mathworks.com/
%       matlabcentral/fileexchange/35094-knee-point), MATLAB Central File 
%       Exchange. Retrieved March 10, 2020.
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% [2] 	Garg, Vaibhav, Ignacio Santamaria, David Ramirez, and Louis Scharf. 
% 	"Subspace Averaging and Order Determination for Source Enumeration." 
% 	IEEE Transactions on Signal Processing (2019).
%
% [3] 	Santamaría, Ignacio, Louis L. Scharf, Chris Peterson, Michael Kirby, 
% 	and J. Francos. "An order fitting rule for optimal subspace averaging." 
% 	In 2016 IEEE Statistical Signal Processing Workshop (SSP), pp. 1-4. 
% 	IEEE, 2016.

% ------------------------------------------------------------------------
% CREATED:      21/05/2019 by Tim Marrinan
%
% LAST EDITED:  23/09/2019 by Tim Marrinan
%
% NOTES: 
%
% 23/09/2019 T. Marrinan
% I am updating the order-selection rule to mirror the paper and cutting
% out the superfluous stuff.  I have archived the code before the changes
% as GMEB_main_deprecated_23_09_2019.m
% ------------------------------------------------------------------------
%clear all; clc;

function GMEB_Exp003( scenario, nRuns )

%% Experimental setup:
%scenario    = 'Fig. 8';                             % Choose scenario
%scenario    = 'Fig. 9';                             % Choose scenario
%scenario    = 'Fig. 10';                            % Choose scenario
%nRuns       = 50;                                   % Number of trials
[param, ~, ~, indVar] =  ...
    GMEB_ScenarioSpecification( scenario );         % Set parameters
k           = param.common_dims_bb;                 % Optimal dimension
kmax        = param.max_dims;                       % Max dimension
nVar        = size(indVar,2);                         % Number of ind. vars.
printout    = 'verbose';                            % 'verbose'/'none'

% Initialize variables
detected_rank.proposed = zeros(nRuns,nVar);
detected_rank.santamaria = zeros(nRuns,nVar);
detected_rank.hybrid = zeros(nRuns,nVar);
detected_rank.svd = zeros(nRuns,nVar);

%% Main loop
for iter = 1 : nRuns
    for ival = 1:nVar
        % Reset the parameters each time
        % (The nonuniform option is a bit wonky and adds more points 
        % each time. I should probably fix that, huh?)
        [param, dual_options, ~, indVar] =  ...
                GMEB_ScenarioSpecification( scenario );    % Set parameters
        switch lower(scenario)            
            case 'fig. 8'
                param.ambient_dims = indVar(ival);
            case 'fig. 9'
                param.noise_var = indVar(ival);
            case 'fig. 10'
                param.ambient_dims = indVar(ival);
            otherwise
                error('Unknown scenario');   
        end

        switch lower(printout)
            case 'verbose'
                clc
                fprintf('\t--------------------------------------------------------\n');
                fprintf('\tExperiment:\t\t\t\t003\n');
                fprintf('\tScenario:\t\t\t\t%s\n',scenario);
                fprintf('\tTrial:\t\t\t\t\t%d/%d\n',iter,nRuns);
                fprintf('\tSNR:\t\t\t\t\t%.4f\n',10*log10(param.common_dims_bb/param.noise_var));
                fprintf('\tAmbient Dimension:\t\t%d\n',param.ambient_dims);
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
        
        % Generate data
        [ data, param ] = GMEB_DataGen( param ); % NestedBall_DataGen_Modified( param );% GMEB_DataGen( param );
        small_k         = param.common_dims_sb;
        
        if param.int_pts
            pts = data.noisy_samples_both;
        else
            pts = data.noisy_samples_all;
        end     
        n_sets = size(pts,1);

        % Run order-selection rule
        dual_options.kmax = kmax;
        [penalty, hybrid, ~] = GMEB_OrderSelection(pts, dual_options);
        
        % Proposed order-selection rule
        detected_rank.proposed(iter,ival) = penalty.kStar;
                
        % Hybrid order-selection rule
        detected_rank.hybrid(iter,ival) = hybrid.kStar;

        % For the Santamaria method [1]
        evals = svd([pts{:}],'econ'); 
        evals = (evals.^2)./n_sets;
        test2 = evals;
        if min(test2)<.5
            evals_dim = find(evals < .5,1);
        else
            evals_dim = size(test2,1)+1;
        end
        detected_rank.santamaria(iter,ival) = evals_dim-1;

        % SVD
        kp = knee_pt(evals);
        detected_rank.svd(iter,ival) = kp - 1;
    end
end

% Compute accuracy and mean selected order
clear accuracy med
accuracy{1} = zeros(nVar,1);
accuracy{2} = zeros(nVar,1);
accuracy{3} = zeros(nVar,1);
accuracy{4} = zeros(nVar,1);
med{1} = zeros(nVar,1);
med{2} = zeros(nVar,1);
med{3} = zeros(nVar,1);
med{4} = zeros(nVar,1);
for i = 1:nVar  
    accuracy{1}(i) = sum(detected_rank.santamaria(1:iter,i) == min(k,small_k))/iter;
    accuracy{2}(i) = sum(detected_rank.proposed(1:iter,i) == min(k,small_k))/iter;
    accuracy{3}(i) =  sum(detected_rank.svd(1:iter,i) == min(k,small_k))/iter;
    accuracy{4}(i) = sum(detected_rank.hybrid(1:iter,i) == min(k,small_k))/iter;
   
    med{1}(i) = mean(detected_rank.santamaria(1:iter,i));
    med{2}(i) = mean(detected_rank.proposed(1:iter,i));
    med{3}(i) = mean(detected_rank.svd(1:iter,i));
    med{4}(i) = mean(detected_rank.hybrid(1:iter,i));
end

%% Performance plots
MethodColors    = cell(6,1);
MethodColors{1} = [251,128,114]/255;
MethodColors{2} = [152,78,163]/255;
MethodColors{3} = [253,180,98]/255;
MethodColors{4} = [141,211,199]/255;

LineStyle{1} = '-o';
LineStyle{2} = '--^';
LineStyle{3} = '-.x';
LineStyle{4} = ':s';

l_width = 3;
fig_width = 12;
fig_height = 12;

f       = cell(3,1);
f{1}    = figure('units','centimeters','Position',...
            [3,3,fig_width,fig_height]);
f{2}    = figure('units','centimeters','Position',...
            [3,3,fig_width,fig_height]);

clear nam
nam{4} = 'SVD';
nam{2} = 'Santamaria et al.';
nam{3} = 'Proposed';
nam{1} = 'k^* = 3';
nam{5} = 'Hybrid';

switch lower(scenario)
    case 'fig. 8'
        % Fig. 8a
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{1})
        for i = 1 : 4
            g{i} = plot(indVar,accuracy{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);
            hold on
        end
        legend({nam{2}, nam{3}, nam{4}, nam{5}})
        figtitle = strcat(scenario,'a');
        title(figtitle)
        ylabel('Accuracy')
        xlabel('Ambient dimension, $n$','interpreter','latex')
        ylim([0,1])
        xlim([20,200])
        hold off

        % Fig. 8b
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{2})
        g{5} = plot(indVar,min([k,small_k]).*ones(nVar,1),'k','LineWidth',3);
        hold on
        for i = 1 : 4
            g{i} = plot(indVar,med{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);    
            hold on
        end
        legend(nam)
        figtitle = strcat(scenario,'b');
        title(figtitle)
        ylabel('Mean selected order')
        xlabel('Ambient dimension, $n$','interpreter','latex')
        ylim([-1,ceil(max(max([med{:}])))])
        xlim([20,200])
        hold off   
        
    case 'fig. 9'
        % Fig. 9a
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{1})
        for i = 1 : 4
            g{i} = plot(10*log10(k./indVar),accuracy{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);
            hold on
        end
        legend({nam{2}, nam{3}, nam{4}, nam{5}})
        figtitle = strcat(scenario,'a');
        title(figtitle)
        ylabel('Accuracy')
        xlabel('SNR (dB)')
        ylim([0,1])
        hold off

        % Fig. 9b
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{2})
        g{5} = plot(10*log10(k./indVar),min([k,small_k]).*ones(nVar,1),'k','LineWidth',3);
        hold on
        for i = 1 : 4
            g{i} = plot(10*log10(k./indVar),med{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);    
            hold on
        end
        legend(nam)
        figtitle = strcat(scenario,'a');
        title(figtitle)
        ylabel('Mean selected order')
        xlabel('SNR (dB)')
        ylim([0,param.max_dims])
        hold off
    case 'fig. 10'
        % Fig. 10a
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{1})
        for i = 1 : 4
            g{i} = plot(indVar,accuracy{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);
            hold on
        end
        legend({nam{2}, nam{3}, nam{4}, nam{5}})
        figtitle = strcat(scenario,'a');
        title(figtitle)
        ylabel('Accuracy')
        xlabel('Ambient dimension, $n$','interpreter','latex')
        ylim([0,1])
        hold off

        % Fig. 10b
        g = cell(size(nam,2),1);
        set(0,'CurrentFigure',f{2})
        g{5} = plot(indVar,min([k,small_k]).*ones(nVar,1),'k','LineWidth',3);
        hold on
        for i = 1 : 4
            g{i} = plot(indVar,med{i},LineStyle{i},'Color',MethodColors{i},'LineWidth',l_width);    
            hold on
        end
        legend(nam)
        figtitle = strcat(scenario,'b');
        title(figtitle)
        ylabel('Mean selected order')
        xlabel('Ambient dimension, $n$','interpreter','latex')
        ylim([-1,ceil(max(max([med{:}])))])
        hold off    
    otherwise
        error('Unknown scenario');   
end