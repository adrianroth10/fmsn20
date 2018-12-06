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
% data
I = ~isnan(Y);
Yvalid = false(5000,1);
for i = 1:10:length(I)
    index = i;
    while I(index)==0 
    index = index+1;
    end
    Yvalid(i)= true;
    I(index) = 0;
end

%create Q-matrix
[u1, u2] = ndgrid(1:sz(1),1:sz(2));
[C,G,G2] = matern_prec_matrices([u1(:) u2(:)]);
%mean value-vector (might not need all)
Bgrid = [ones(prod(sz),1) bei_grad(:).^2 bei_elev(:) bei_elev(:).^2 ];
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
G, G2, qbeta, true), [0 0]);
%conditional mean is given by the mode
E_xy = x_mode;
E_zy = exp(Atilde * E_xy);

%reuse taylor expansion to compute posterior precision
tau = par(1);
kappa2 = par(2);
Qcar = tau*(kappa2*C + G);
Qsar = tau*(kappa2^2*C + 2*kappa2*G + G2);
Nbeta = size(Bgrid,2);
Qtilde = blkdiag(Qcar, qbeta*speye(Nbeta));

[~, f, Q_xy] = GMRF_taylor(E_xy, Y(I), Atilde(I, :), Qtilde);

%% Variance of beta
e = [zeros(size(Q_xy,1)-size(Bgrid,2), size(Bgrid,2)); eye(size(Bgrid,2))];
V_beta0 = e'*(Q_xy\e);


%1000 samples from the approximate posterior
Rxy = chol(Q_xy);
x_samp = bsxfun(@plus, E_xy, Rxy\randn(size(Rxy,1),1000));

Vzy = zeros(5000,1);
for i =1:5000
    Vzy(i) = var(x_samp(i,:));
end
%Vzy = Atilde*V_beta0*Atilde';

rms_error = sqrt(mean(~isnan(Vzy(Yvalid).*(E_zy(Yvalid)- Y(Yvalid)).^2)));

imagesc(reshape(Vzy,sz))
%% Plotting
figure()
subplot(2,2,1)
title('Counted data')
imagesc(reshape(bei_counts, sz))
colorbar
subplot(2,2,2)
title('Estimated data')
imagesc(reshape(E_zy, sz))
%hold on 
%scatter(E_zy(Yvalid), E_zy(Yvalid), 20, sqrt(Vzy(Yvalid)), 'filled','markeredgecolor','g')
colorbar
subplot(2,2,3)
title('Elevation')
imagesc(reshape(bei_elev,sz))
colorbar
subplot(2,2,4)
title('Elevation gradient')
imagesc(reshape(bei_grad,sz))
colorbar
