function Y = GrExp(X,TY)
% SYNTAX:   [Y] = GrExp(X,TY)
%
% INPUTS:   'X' is the point about which the tangent space has been
%           computed.
%
%           'TY' is a point in the tangent space of X.
%
% OUTPUTS:  'Y' is the representation of the point TY on the Grassmannian
%           manifold.
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


[U, S, V] = svd(TY,0);
Y = X*V*diag(cos(diag(S))) + U*diag(sin(diag((S))));
