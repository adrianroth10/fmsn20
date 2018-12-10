%load data
load HA2_forest
%size of the grid
sz = size(bei_counts);
%observations
Y = bei_counts(:);
% data
rng(0)  % setting seed for predictable random sequence
I = ~isnan(Y);
ind_i = find(I);
Ivalid = logical(zeros(size(I)));
Ivalid(ind_i(rand(size(ind_i)) < 0.1)) = true;
I(Ivalid) = false;

%create Q-matrix
[u1, u2] = ndgrid(1:sz(1),1:sz(2));
[C,G,G2] = matern_prec_matrices([u1(:) u2(:)]);
%mean value-vector (might not need all)
%and observation matrix for the grid
Agrid = speye(prod(sz));
%Agrid = 0*speye(prod(sz));
% Bgrid = [ones(prod(sz),1) bei_elev(:)];
% Bgrid = [ones(prod(sz),1) bei_grad(:)];
% % Bgrid = [ones(prod(sz),1) bei_elev(:) bei_grad(:)];
% Bgrid = [ones(prod(sz),1) bei_elev(:) bei_elev(:).^2];
% Bgrid = [ones(prod(sz),1) bei_grad(:) bei_grad(:).^2];
Bgrid = [ones(prod(sz),1) bei_elev(:) bei_elev(:).^2 bei_grad(:) bei_grad(:).^2];
Nbeta = size(Bgrid,2);
qbeta = 1e-6;
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
G, G2, qbeta, false), [0 0]);
%conditional mean is given by the mode
E_xy = x_mode;
E_zy = exp(Atilde * E_xy);

%reuse taylor expansion to compute posterior precision
tau = exp(par(1));
kappa2 = exp(par(2));
gamma = 1;
Qcar = tau*(kappa2*C + G);
Qsar = tau*(kappa2^2*C + 2*kappa2*G + G2);
Qosc = tau*(kappa2^2*C + 2*gamma*kappa2*G + G2);
Nbeta = size(Bgrid,2);
Qtilde = blkdiag(Qsar, qbeta*speye(Nbeta));

[~, f, Q_xy] = GMRF_taylor(E_xy, Y(I), Atilde(I, :), Qtilde);

%% Variance of beta
e = [zeros(size(Q_xy,1)-size(Bgrid,2), size(Bgrid,2)); eye(size(Bgrid,2))];
V_beta0 = e'*(Q_xy\e);
beta = E_xy(end - length(V_beta0) + 1:end)
beta_int_size = sqrt(diag(V_beta0)) * 1.96


%1000 samples from the approximate posterior
Rxy = chol(Q_xy);
x_samp = bsxfun(@plus, E_xy, Rxy\randn(size(Rxy,1),1000));


Vx  = 1/(size(x_samp,2)-length(V_beta0))*sum(x_samp(1:end-length(V_beta0),:).^2,2);
Vzy = Vx + sum((Bgrid*V_beta0).*Bgrid,2);
std = Vzy.^(1/2);

rms_error = validations(Y(Ivalid), E_zy(Ivalid), Vzy(Ivalid));  % sqrt(mean(Vzy(Ivalid).*(E_zy(Ivalid)- Y(Ivalid)).^2)));

type = 1; % car model
proj2_run;

%%
% figure()
% subplot(2,2,1)
% title('Counted data')
% imagesc(reshape(bei_counts, sz))
% colorbar
% subplot(2,2,2)
% title('Counted data')
% imagesc(reshape(Ivalid, sz))
% colorbar
% subplot(2,2,3)
% title('Elevation')
% imagesc(reshape(bei_elev,sz))
% colorbar
% subplot(2,2,4)
% title('Elevation gradient')
% imagesc(reshape(bei_grad,sz))
% colorbar
