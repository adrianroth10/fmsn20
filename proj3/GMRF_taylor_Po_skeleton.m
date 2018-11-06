function [f, df, d2f]= GMRF_taylor_Po(x_0, y, A, Q, mu)
% GMRF_TAYLOR_PO_SKELETON  Taylor expansion of the conditional for non-Gaussian observations
%
% [f, df, d2f]= GMRF_taylor_Po(x_0, y, A, Q, mu)
%
% x_0 = value at which to compute taylor expansion
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
% mu = the mean value vector
%
% Function should return taylor expansion, gradient and Hessian.
%
% This is only a skeleton for Home Assignment 3.

% $Id: gmrf_negloglike_skeleton.m 4454 2011-10-02 18:29:12Z johanl $

% Remove this line from your copy:
%warning('This is only a skeleton function!  Copy it and fill in the blanks!')

%compute the function
logp = [];
f = [];

%compute derivatives
d_logp = [];
df = [];

%compute hessian
d2_logp = [];
n = size(Q,1);
d2f = Q - spdiags(A'*d2_logp,0,n,n);

