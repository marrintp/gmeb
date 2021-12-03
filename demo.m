%   File: demo.m
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
% This is a script not a function.  Parameters are set within. 
%
% ------------------------------------------------------------------------
% OVERVIEW:
% This script generates subspaces of differing dimesions and finds the GMEB
% center of optimal dimension. It then embeds that data and average in 2
% or 3 dimensions and displays the embedding.
%
% ------------------------------------------------------------------------
% INPUTS:
% n/a
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
% [c] GMEB_OrderSelection.m
%       Not included. 
%       Written by Tim Marrinan, see file for documentation.
%
% [d] GMEB_DualSubgrad.m       
%       Not included.
%       Written by Tim Marrinan, see file for documentation.
%
% [e] diffusionKernel.m       
%       Not included.
%       See file for documentation.
%
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% ------------------------------------------------------------------------
% CREATED:      21/03/2020 by Tim Marrinan
%
% LAST EDITED:  21/03/2020 by Tim Marrinan
%
% NOTES: 
% 21/03/2020
% The embedding is never going to actually look like a ball unless the
% points are all 1-dimensional subspaces on a sphere (tm)
%
% ------------------------------------------------------------------------
path1 = genpath('src');
path2 = genpath('examples');
addpath(path1);
addpath(path2);

scenario    = 'custom';                             % Choose scenario
[param, dual_options, ~, ~] =  ...
    GMEB_ScenarioSpecification( scenario );         % Set parameters

% Parameters for visualizing the output
method      = 'multidimensional scaling';           % 'multidimensional scaling'/'diffusion map'/'pca'
k           = param.common_dims_bb;                 % Optimal dimension
d_type      = param.d_type;                         % Distance measure
printout    = 'verbose';                            % 'verbose'/'none'
threeD      = false;                                % true/false;



%% Main script
switch lower(printout)
    case 'verbose'
        clc        
        fprintf('\n\tGMEB demo\n');
        fprintf('\t--------------------------------------------------------\n');
        fprintf('\tScenario:\t\t\t\t%s\n',scenario);
        fprintf('\tDim. reduction:\t\t\t%s\n',method);
        fprintf('\tAmbient dimensions:\t\t%g\n',param.ambient_dims);
        fprintf('\tSNR:\t\t\t\t\t%g\n',10*log10(k./param.noise_var));
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
[ data, param ] = GMEB_DataGen( param );


% Specify point set from generated data
if param.int_pts && param.nonuni
    pts = data.noisy_samples_both;
    nBB = ceil(param.num_samples_bb*(param.int_pts_ratio+1));
    nSB = ceil(param.num_samples_sb*(param.int_pts_ratio+1));
elseif param.int_pts && ~param.nonuni
    pts = data.noisy_samples_both;
    nBB = ceil(param.num_samples_bb*(param.int_pts_ratio+1));
    nSB = ceil(param.num_samples_sb*(param.int_pts_ratio+1));    
elseif ~param.int_pts && param.nonuni
    pts = data.noisy_samples_all;
    nBB = param.num_samples_bb;
    nSB = param.num_samples_sb;
elseif ~param.int_pts && ~param.nonuni
    pts = data.noisy_samples_all;
    nBB = param.num_samples_bb;
    nSB = param.num_samples_sb;
end
n_sets = size(pts,1);


%% Select optimal order
[penalty, hybrid, dual_options] = GMEB_OrderSelection(pts, dual_options);

switch lower(printout)
    case 'verbose'
        fprintf('\tGround truth order:\t\t%g\n',k);
        fprintf('\tPenalty order:\t\t\t%g\n',penalty.kStar);
        fprintf('\tHybrid order:\t\t\t%g\n',hybrid.kStar);
        fprintf('\t--------------------------------------------------------\n');
    otherwise
end

%% Compute the minimum enclosing ball 
[iterates, final, dual_options] = GMEB_DualSubgrad(pts, k, dual_options);
GMEB = final.u;


%% Embed data into R^2 or R^3
pts{n_sets+1} = GMEB;
pts{n_sets+2} = data.center;               
distMat = zeros(size(pts,1),size(pts,1));
for i = 1 : size(pts,1)
    for j = i + 1 : size(pts,1)
        distMat(i,j) = GrDist(pts{i},pts{j},k,d_type);
        distMat(j,i) = distMat(i,j);
    end
end
switch lower(method)
    case 'multidimensional scaling'
        clear coords evals L_sharp
        [coords,~] = cmdscale(distMat,3);
        cc_GMEB = coords(n_sets+1,:);
        cc_true = coords(n_sets+2,:);
    case ' diffusion map'
        clear coords evals L_sharp
        [coords] = diffusionKernel(distMat,20,1,3);
        cc_GMEB = coords(n_sets+1,:);
        cc_true = coords(n_sets+2,:);
    case 'pca'
        % Need to make everything point the same way
        [~,ind] = max(abs(pts{n_sets+3}));
        for i = 1 : n_sets+3
            if pts{i}(ind) < 0
                j = j + 1;
                pts{i} = -1*pts{i};
            end
        end
        coeff = pca([pts{:}]');
        coords = [pts{:}]'*coeff;
        cc_GMEB = coords(n_sets+1,:);
        cc_true = coords(n_sets+2,:);
    otherwise
        error('Unknown embedding method');   
end

    
%% Visualize embedding   
true_center_color = [17,119,51]/255;
gmeb_color = [221 204 119]/255;

fig_width = 12;
fig_height = 12;
figure('units','centimeters','Position',[4,4,fig_width,fig_height]);


if threeD
    switch lower(printout)
        case 'verbose'
            fprintf('\n\tNote:\tAntipodal points on the sphere represent\n');
            fprintf('\t\t\tthe same point on Gr(1,3), so only one hemisphere\n');
            fprintf('\t\t\twill appear to contain points...\n');
        otherwise
    end
    exes = scatter3(coords(1:nBB,1),coords(1:nBB,2),coords(1:nBB,3),...
        'x','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],...
        'DisplayName','$\*X_i \in \mathcal{B}_{0.5}(\mathbf{Z}_1)$');
    hold on
    ohs = scatter3(coords(nBB+1:nBB+nSB,1),coords(nBB+1:nBB+nSB,2),...
        coords(nBB+1:nBB+nSB,3),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 1 1],'DisplayName',...
        '$\*X_i \in \mathcal{B}_{0.05}(\mathbf{Z}_2)$');
    tm = scatter3(cc_true(1)',cc_true(2)',cc_true(3)',200,'s',...
        'MarkerEdgeColor','k','MarkerFaceColor',true_center_color,...
        'DisplayName','True center');
    gm = scatter3(cc_GMEB(1)',cc_GMEB(2)',cc_GMEB(3)',100,'d',...
        'MarkerEdgeColor','k','MarkerFaceColor',gmeb_color,...
        'DisplayName','Estimated center');
    title(scenario)
    xlim([min(coords(:,1)),max(coords(:,1))])
    ylim([min(coords(:,2)),max(coords(:,2))])
    zlim([min(coords(:,3)),max(coords(:,3))])
    l = legend('Location','best');
    set(l, 'Interpreter', 'latex')
    set(gca,'visible','off')
    hold off
    %alpha(km,.75)
    %alpha(gm,.75)
else
    exes = scatter(coords(1:nBB,1),coords(1:nBB,2),'x',...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],...
        'DisplayName','$\*X_i \in \mathcal{B}_{0.5}(\mathbf{Z}_1)$');
    hold on
    ohs = scatter(coords(nBB+1:nBB+nSB,1),coords(nBB+1:nBB+nSB,2),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'DisplayName',...
        '$\*X_i \in \mathcal{B}_{0.05}(\mathbf{Z}_2)$');
    tm = scatter(cc_true(1)',cc_true(2)',200,'s',...   
        'MarkerEdgeColor','k','MarkerFaceColor',true_center_color,...
        'DisplayName','True center');
    gm = scatter(cc_GMEB(1)',cc_GMEB(2)',100,'d',...
        'MarkerEdgeColor','k','MarkerFaceColor',...
        gmeb_color,'DisplayName','Estimated center');
    title(scenario)
    xlim([min(coords(:,1)),max(coords(:,1))])
    ylim([min(coords(:,2)),max(coords(:,2))])
    l = legend('Location','best');
    set(l, 'Interpreter', 'latex')
    set(gca,'visible','off')
    hold off
end






