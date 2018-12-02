%load data
load HA2_forest

subplot(221)
imagesc( bei_counts )
subplot(222)
imagesc( isnan(bei_counts) )
subplot(223)
imagesc( bei_elev )
subplot(224)
imagesc( bei_grad )

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
figure
subplot(121)
spy(G2)
subplot(122)
spy(G2(p,p))

%reorder precision matrices
C = C(p,p);
G = G(p,p);
G2 = G2(p,p);
%and observation matrix
Agrid = Agrid(:,p);

%create A tilde matrix
Atilde = [Agrid Bgrid]; 
%%
%we need a global variable for x_mode to reuse
%between optimisation calls
global x_mode;
x_mode = [];

%subset Y and Atilde to observed points
par = fminsearch( @(x) GMRF_negloglike(x, Y(I), Atilde(I,:), C, G, G2, 1e-6, true), [0 0]);
%conditional mean is given by the mode  
E_xy = x_mode;
%and reconstruction (field+covariates)
E_zy = Atilde*x_mode;
imagesc( reshape(E_zy,sz) )
%%
%reuse taylor expansion to compute posterior precision
[~, ~, Q_xy] = GMRF_taylor(E_xy, Y(I),  Atilde(I,:), Qtilde);

%1000 samples from the approximate posterior
Rxy = chol(Q_xy);
x_samp = bsxfun(@plus, E_xy, Rxy\randn(size(Rxy,1),1000));