function negloglike = GMRF_negloglike_Po(theta, y, Aall, C, G, G2, Qbeta, p)
% GMRF_NEGLOGLIKE_PO_SKELETON  Calculate the GMRF data likelihood, non-Gaussian observations
%
% negloglike = GMRF_negloglike_Po(theta, y, A, C, G, G2, B)
%
% theta = [log([tau2; kappa;]); beta]
% y = the data vector, as a column with n elements
% Aall = the observation matrix, sparse n-by-(N+Nbeta)
% C,G,G2 = matrices used to build a Matérn-like precision,
%          see matern_prec_matrices, sparse N-by-N
% Qbeta = Precision matrix for the regression parameters (N+Nbeta)-by-(N+Nbeta)
% p = reordering vector for the full Q matrix, (N+Nbeta)-by-1
%
% This is only a skeleton for Home Assignment 3.

% $Id: GMRF_negloglike_Po_skeleton.m 4795 2014-11-27 16:45:06Z johanl $

% Remove this line from your copy:
%error('This is only a skeleton function!  Copy it and fill in the blanks!')

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

%set optimisation options
%(gradient, hessian, and reduced accuracy)
if verLessThan('matlab', '8.1')
  %optimoptions was introduces in matlab 2013a (8.1)
  opts = optimset(optimset('fminunc'), 'Largescale', 'on', ...
                  'GradObj', 'on', 'Hessian', 'on', ...
                  'Display', 'off', 'TolFun', 1e-3, 'TolX', 1e-3);
else
  %modern matlab, use latest functions
  opts = optimoptions('fminunc', 'Algorithm', 'trust-region', ...
                      'GradObj', 'on', 'Hessian', 'on', ...
                      'Display', 'off', 'TolFun', 1e-3, 'TolX', 1e-3);
end
%find mode
x_mode = fminunc(@(x) GMRF_taylor_Po(x, y, Aall, Qall), x_mode, opts);

%find the Laplace approximation of the denominator
[logp, ~, Q_xy] = GMRF_taylor_Po(x_mode, y, Aall, Qall);
%note that logp = -log_obs + x_mode'*Q*x_mode/2.

%Compute choleskey factors
[R_x,p_x] = chol(Q_x);
[R_xy,p_xy] = chol(Q_xy);
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
