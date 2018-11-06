function [logp, D_logp, D2_logp]= GMRF_taylor_Po(x_0, y, A, Q)
% GMRF_TAYLOR_PO_SKELETON  Taylor expansion of the conditional for non-Gaussian observations
%
% [logp, D_logp, D2_logp]= GMRF_taylor_Po(x_0, y, A, Q)
%
% x_0 = value at which to compute taylor expansion
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
%
% Function should return taylor expansion, gradient and Hessian.
%
% This is only a skeleton for Home Assignment 2.

% $Id: GMRF_taylor_Po_skeleton.m 4792 2014-11-27 14:30:32Z johanl $

% Remove this line from your copy:
error('This is only a skeleton function!  Copy it and fill in the blanks!')

%compute log observations p(y|x)
z = A*x_0;
f = [];

%compute log p(x|y,theta)
logp = [];

if nargout>1
  %compute derivatives (if needed, i.e. nargout>1)
  df = [];
  D_logp = [] - A'*df;
end

if nargout>2
  %compute hessian (if needed, i.e. nargout>2)
  d2f = [];
  n = size(A,1);
  D2_logp = [] - A'*spdiags(d2f,0,n,n)*A;
end
