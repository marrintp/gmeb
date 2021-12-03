function [MU, err, iters, tellapsed] = KarcherMean(data, start_pt, tol, maxiters,check)
% SYNTAX:   [MU, err, iters, tellapsed] = KarcherMean(data, start_pt, tol, maxiters,check)
%
% INPUTS:   'data' is the collection of subspaces you wish to find a
%           representative for.  It should be entered as a num_pts by 1
%           cell array with each cell being an n by q_i matrix, where q_i 
%           CANNOT vary.  In other words, the datacubes need to be unrolled
%           already when using this function
%
%           'start_pt' is a point from 'data' that you would like to use as
%           your initial base point.
%
%           'tol' how close MU needs to be at step n, relative to MU at
%           step n+1 to say that the algorithm has converged.
%
%           'maxiters' is the maximum allowable number of iterations.  The
%           algorithm will stop iterating if this threshold is reached and
%           return the current value of MU.
%
%           'check' determines whether or not the algorithm will verify
%           that your data points are in fact orthonormal matrices.  If
%           check == 1, it does the check, otherwise it does not.
%
% OUTPUTS:  'MU' is the Karcher mean, returned as a matrix. 
%
%           'err' is the distance between MU at the second to last 
%           iteration and MU at the last iteration.
%
%           'tellapsed' is the total time it took to compute MU.
%
% NOTES:    This function calls GrExp, GrLog, and GrDist which are NOT
%           included.
%
% LAST EDITED: 04/28/2013 by Tim Marrinan
%
%--------------------------------------------------------------------------
% REFERENCE:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, J. Ross Beveridge, Bruce Draper, Michael Kirby, and 
%	Chris Peterson. "Finding the subspace mean or median to fit your 
% 	need." In Proceedings of the IEEE Conference on Computer Vision and 
%	Pattern Recognition, pp. 1082-1089. 2014.
%
%--------------------------------------------------------------------------

MU = 0;
tellapsed = 0;
etol = 10e-10;
err = 999999;
iters = 1;

tic;
tstart = tic;

%test to make sure points are on Grassmannian
[nr nc] = size(start_pt);
NumAngles = nc;
M = start_pt'*start_pt;
II = eye(nc);
if norm(M-II,'fro') < etol
else
    display('Starting point is not orthonormal.');
    [U S V] = svd(start_pt,0);
    start_pt = U*V';
end

[num_pts other] = size(data);
[~, ncz] = size(data{1});
if other>1
    display('Halting due to dimension mismatch.');
    return;
end
if check == 1
    for i = 1:num_pts
        Y = squeeze(data{i});
        M = Y'*Y;
        II = eye(ncz);
    
        if norm(M-II,'fro') < etol
        %display('Z already on grassmannian');
        else%find closest o.n. matrix
        %display('Finding closest on matrix')
        [U, ~, V] = svd(Y,0);
        data{i} = U*V';
        end
    end
end

% Begin main loop of the Karcher mean
MU = start_pt;
while err > tol
    %iters
    if mod(iters,100) == 0
        %iters
        %err
    end
    SUMT = zeros(nr, nc);
    
    for i = 1:num_pts
       SUMT  = SUMT + GrLog(MU,data{i});
    end
    SUMT = SUMT/num_pts;
    
    NewMU = GrExp(MU, SUMT);
    olderr = err;
    [err] = GrDist(NewMU, MU, size(MU,2), 'DGEO');
    err;

    MU = NewMU;
    
    if err < tol
        %display('Converged')
        tellapsed = toc(tstart);
        break;
    end
    iters = iters  + 1;
    if iters >= maxiters
        tellapsed = toc(tstart);
        return
    end
end



