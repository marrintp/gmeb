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
%           I commented out the code that finds an orthonormal basis for
%           each point because I realized that the qr-factorization makes
%           Q n by n even if the input is not rank n.
%
% LAST EDITED: 03/13/13 by Tim Marrinan
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

% Set up the default values
distance = 500;

if ~exist('dist_type','var')
    dist_type = 'DGEO';
end

q = min(size(Q,2),size(P,2));
offset = abs(size(Q,2)-size(P,2));
N=min(N,q);

%Euclidean
if strcmp(dist_type,'EUC')
    distance = norm(P-Q,'fro');
    return
end

% % Test to make sure the matrices are orthonormal
dd=norm(P'*P-eye(size(P,2)));
ee=norm(Q'*Q-eye(size(Q,2)));

if dd>100*eps
    %disp('Finding nearest pt. on manifold.')
    [P, ~]=qr(P,0);
end

if ee>100*eps
    %disp('Finding nearest pt. on manifold.')
    [Q, ~]=qr(Q,0);
end

if size(P,1)~=size(Q,1)
    disp('GrDist halting due to dimension mismatch.')
    size(P,1)
    size(Q,1)
    return
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



