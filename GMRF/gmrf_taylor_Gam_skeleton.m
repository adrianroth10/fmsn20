function [logp, Dlogp, D2logp]= gmrf_taylor_Gam(x_0, y, A, Q, b)
% GMRF_TAYLOR_GAM Taylor expansion of the conditional for non-Gaussian observations
%
% [f, df, d2f]= GMRF_taylor_Gam(x_0, y, A, Q, b)
%
% x_0 = value at which to compute taylor expansion
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
% b = variance parameter of Gamma distribution
%
% Function should return taylor expansion, gradient and Hessian.
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_skeleton.m 4454 2011-10-02 18:29:12Z johanl $
error('This is only a skeleton function!  Copy it and fill in the blanks!')

%Observation model
% Y_i ~ Gamma(b, exp(z_i)/b)

%compute log observations, and derivatives
z = A*x_0;
%compute log observations p(y|x)
%this should include all constants containing b!
f = [];

%compute -log p(x|y,theta) (negative since we're doing minimisation)
logp = [];

if nargout>1
  %compute derivatives (if needed, i.e. nargout>1)
  df = [];
  Dlogp = [] - A'*df;
end

if nargout>2
  %compute hessian (if needed, i.e. nargout>2)
  d2f = []
  n = size(A,1);
  D2logp = [] - A'*spdiags(d2f,0,n,n)*A;
end
