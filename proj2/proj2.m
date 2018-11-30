%%1. Start by computing all the components in the theory section
%%2. Complete gmrf taylor skeleton.m.
%%3. Complete gmrf negloglike skeleton.m.
%%4. Estimate parameters for the different models and decide on covariates.
%%5. Reconstruct the different components of the latent field, e.g. E(z|y),
%%and compute the reconstruction uncertainty, V(z|y).

%% Suggested skeleton
%load data
load HA2_forest
%size of the grid
sz = size(bei_counts);
%observations
Y = bei_counts(:);
%missing data
I = ~isnan(Y);

%create Q-matrix
[u1, u2] = ndgrid(1:sz(1),1:sz(2));
[C,G,G2] = matern_prec_matrices([u1(:) u2(:)]);
%mean value-vector (might not need all)
Bgrid = [ones(prod(sz),1) bei_elev(:) bei_grad(:)];
%and observation matrix for the grid
Agrid = speye(prod(sz));
%G2 is the most dense of the matrices, lets reorder
p = amd(G2);
%reorder precision matrices
C = C(p,p);
G = G(p,p);
G2 = G2(p,p);
%and observation matrix
Agrid = Agrid(:,p);
%create A tilde matrix
Atilde = [Agrid Bgrid];
%we need a global variable for x_mode to reuse
%between optimisation calls
global x_mode;
x_mode = [];
%subset Y and Atilde to observed points
par = fminsearch( @(x) GMRF_negloglike(x, Y(I), Atilde(I,:), C, ...
G, G2, 1e-6, true), [0 0]);
%conditional mean is given by the mode
E_xy = x_mode;
%reuse taylor expansion to compute posterior precision
tau = 1;
kappa2 = 0.001;
Qcar = tau*(kappa2*C + G);
Qsar = tau*(kappa2^2*C + 2*kappa2*G + G2);
[~, ~, Q_xy_SAR] = GMRF_taylor(E_xy, Y(I), Atilde, Qcar);


e = [zeros(size(Q_xy,1)-size(B,2), size(B,2)); eye(size(B,2))];
V_beta0 = e'*(Q_xy\e);


