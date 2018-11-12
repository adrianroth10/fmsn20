function [logp, D_logp, D2_logp]= GMRF_taylor_Po(x_0, y, A, Q, E)
% GMRF_TAYLOR_PO_SKELETON  Taylor expansion of the conditional for non-Gaussian observations
%
% [logp, D_logp, D2_logp]= GMRF_taylor_Po(x_0, y, A, Q)
%
% x_0 = value at which to compute taylor expansion, as a column with N elements
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
% E = the population count in each region, as a column with n elements
%
% Function should return taylor expansion, gradient and Hessian.
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_taylor_Po_skeleton.m 4897 2015-11-22 11:36:48Z johanl $

% Remove this line from your copy:
error('This is only a skeleton function!  Copy it and fill in the blanks!')

%ensure that E=1 if E not given (i.e. same/no population weight in all regions)
if nargin<5, E=1; end

%compute log observations p(y|x)
z = A*x_0;
f = [];

%compute -log p(x|y,theta)
logp = [] - sum(f);

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
