function [TY] = GrLog(X,Y)
% SYNTAX:   [TY] = GrLog(X,Y)
%
% INPUTS:   'X' is the point about which the tangent space has been
%           computed.
%
%           'Y' is the point on the Grassmannian manifold that is to be
%           mapped to the tangent space of X.
%
% OUTPUTS:  'TY' is the representation of Y in the tangent space of X.
%
% NOTES:    
%
% LAST EDITED: 02/24/13 by Tim Marrinan
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

m = size(X,1);

% temp = (eye(m)-X*X')*Y*inv(X'*Y);
% The following line is a slightly faster way to compute temp.
% temp = (eye(m)*Y)*inv(X'*Y) - (X*(X'*Y))*inv(X'*Y);
temp = (eye(m)*Y)/(X'*Y) - (X*(X'*Y))/(X'*Y);
[U,S,V] = svd(temp,0);
Theta = atan(S);

TY = U*Theta*V';