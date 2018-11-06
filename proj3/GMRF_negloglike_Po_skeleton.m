function negloglike = GMRF_negloglike_Po(theta, y, A, C, G, G2, B)
% GMRF_NEGLOGLIKE_PO_SKELETON  Calculate the GMRF data likelihood, non-Gaussian observations
%
% negloglike = GMRF_negloglike_Po(theta, y, A, C, G, G2, B)
%
% theta = [log([tau2; kappa;]); beta]
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% C,G,G2 = matrices used to build a Matérn-like precision,
%          see matern_prec_matrices, sparse N-by-N
% B = the expectation basis matrix, N-by-length(beta)
%
% This is only a skeleton for Home Assignment 3.

% $Id: gmrf_negloglike_skeleton.m 4454 2011-10-02 18:29:12Z johanl $

% Remove this line from your copy:
%warning('This is only a skeleton function!  Copy it and fill in the blanks!')

tau2 = exp(theta(1));
kappa = exp(theta(2));
beta = theta(3:end);

%first create Q and mean value
mu = [];
Q_x = [];

%find mode
opts = optimset('GradObj','on', 'Hessian', 'on', 'Display', 'off');
x_mode = fminunc(@(x) GMRF_taylor_Po(x, y, A, Q_x, mu), mu, opts);
%find the Laplace approximation of the denominator
[f, df, Q_xy] = GMRF_taylor_Po(x_mode, y, A, Q_x, mu);

%amd-reorder for faster choleskey
p = amd(Q_x);
Q_x = Q_x(p,p);
Q_xy = Q_xy(p,p);
A = A(:,p);
x_mode = x_mode(p);
mu = mu(p);

%Compute choleskey factor
[R_x,p_x] = chol(Q_x);
[R_xy,p_xy] = chol(Q_xy);
if p_x~=0 || p_xy~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  negloglike = realmax;
  return;
end

negloglike = [];
