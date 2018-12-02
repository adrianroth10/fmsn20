function [g, dg, d2g]= GMRF_taylor(x_0, y, A, Q)
% GMRF_TAYLOR  Taylor expansion of the conditional for non-Gaussian observations
% x_0 = value at which to compute taylor expansion
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
%
% Function should return taylor expansion of
%   g = -log p(y|x) + 1/2 x'*Q*x = -f + 1/2 x'*Q*x
% as well as gradient and Hessian. Sign shift since matlab does
% minimisation instead of maximisation

%compute log observations, and derivatives
z = A*x_0;
% size(A)
% size(y)
%  size(Q)
%  size(x_0)
% size(z)
logp = y.* z - exp(z) - log(factorial(y));
% size(logp)
%compute the function
g = x_0'*Q*x_0/2 - sum(logp);

if nargout>1
  %compute derivatives (if needed, i.e. nargout>1)
  d_logp = y-exp(z);
  dg = Q*x_0 - A'*d_logp;
end

if nargout>2
  %compute hessian (if needed, i.e. nargout>2)
  d2_logp = -exp(z);
  n = size(A,1);
  d2g = Q - A'*spdiags(d2_logp,0,n,n)*A;
end