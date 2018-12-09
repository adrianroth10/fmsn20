function negloglike = GMRF_negloglike(theta, y, Atilde, C, G, G2, qbeta, type)
% GMRF_NEGLOGLIKE  Calculate the GMRF data likelihood, non-Gaussian observations
% theta = [log(tau2) log(kappa2)]
% y = the data vector, as a column with n elements
% Atilde = the observation matrix, sparse n-by-(N+Nbeta)
% C,G,G2 = matrices used to build a Matern-like precision,
%          see matern_prec_matrices, sparse N-by-N
% qbeta = prior precision the regression parameters scalar
% type = 1 (CAR), 2 (SAR), 3 (OSC)

tau = exp(theta(1));
kappa2 = exp(theta(2));

%comput Q for a CAR(1), SAR(1) or OSC SAR process
if type == 1
  Q_x = tau*(kappa2*C + G);
elseif type == 2
  Q_x = tau*(kappa2^2*C + 2*kappa2*G + G2);
else
  gamma = (exp(theta(3)) - 1) / (exp(theta(3)) + 1);
  Q_x = tau*(kappa2^2*C + 2*gamma*kappa2*G + G2);
end

%combine Q_x and qbeta
Nbeta = size(Atilde,2) - size(Q_x,1);
Qtilde = blkdiag(Q_x, qbeta*speye(Nbeta));

%declare x_mode as global so that we start subsequent optimisations from
%the previous mode (speeds up nested optimisation).
global x_mode;
if isempty(x_mode)
  %no existing mode, compute a rough initial guess assuming Gaussian errors of log(y+0.1)
  x_mode = (Qtilde + Atilde'*Atilde)\(Atilde'*log(y+.1));
end

%find mode using Newton-Raphson
x_mode = fminNR(@(x) GMRF_taylor(x, y, Atilde, Qtilde), x_mode);

%find the Laplace approximation of the denominator
[g, ~, Q_xy] = GMRF_taylor(x_mode, y, Atilde, Qtilde);
%HINT: What does g contain?
% g = -log p(y|x) + 1/2 x'*Q*x = -f + 1/2 x'*Q*x =   log(p(x˜|y,?))

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

negloglike = g - sum(log(diag(R_x))) + Nbeta / 2 * log(qbeta) + sum(log(diag(R_xy)));

%print diagnostic/debug information (optimization progress)
% fprintf(1, 'Theta: %11.4e %11.4e; fval: %11.4e\n', theta(1), theta(2), negloglike);
