function negloglike = GMRF_negloglike_Po(theta, y, Aall, C, G, G2, Qbeta, p, E)
% GMRF_NEGLOGLIKE_PO_SKELETON  Calculate the GMRF data likelihood, non-Gaussian observations
%
% negloglike = GMRF_negloglike_Po(theta, y, Aall, C, G, G2, Qbeta, p, E)
%
% theta = log([tau2; kappa])
% y = the data vector, as a column with n elements
% Aall = the observation matrix, sparse n-by-(N+Nbeta)
% C,G,G2 = matrices used to build a Matérn-like CAR/SAR precision,
%          see matern_prec_matrices, sparse N-by-N
% Qbeta = Precision matrix for the regression parameters (N+Nbeta)-by-(N+Nbeta)
% p = reordering vector for the full Q matrix, (N+Nbeta)-by-1
% E = the population count in each region, as a column with n elements
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_Po_skeleton.m 4897 2015-11-22 11:36:48Z johanl $

% Remove this line from your copy:
error('This is only a skeleton function!  Copy it and fill in the blanks!')

%ensure that E=1 if E not given (i.e. same/no population weight in all regions)
if nargin<9, E=1; end

%extract parameters
tau = exp(theta(1));
kappa2 = exp(theta(2));

%comput Q for a CAR(1) or SAR(1) process
Q_x = [];

%combine this Q and Qbeta
Qall = blkdiag(Q_x, Qbeta);
%reorder Qall (Aall should be reorderd befoer calling the function)
Qall = Qall(p,p);

%declare x_mode as global so that we start subsequent optimisations from
%the previous mode (speeds up nested optimisation).
global x_mode;
if isempty(x_mode)
  %no existing mode, compute a rough initial guess assuming Gaussian errors
  x_mode = (Qall + Aall'*Aall)\(Aall'*log(y+.1));
end

%find mode - using Newton-Raphson optimisation
x_mode = fminNR(@(x) GMRF_taylor_Po(x, y, Aall, Qall, E), x_mode);

%find the Laplace approximation of the denominator
[logp, ~, Q_xy] = GMRF_taylor_Po(x_mode, y, Aall, Qall, E);
%note that logp = -log_obs + x_mode'*Q*x_mode/2.

%Compute choleskey factors
[R_x, p_x] = chol(Q_x);
[R_xy, p_xy] = chol(Q_xy);
if p_x~=0 || p_xy~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  negloglike = realmax;
  return;
end

%note that logp = -log_obs + x_mode'*Q*x_mode/2.
negloglike = [];

%print diagnostic information (progress)
fprintf(1, 'Theta: %11.4e %11.4e; fval: %11.4e\n', ...
  theta(1), theta(2), negloglike);
