%% 1.

% dist = 0:80;
% sigma2 = 1;
% kappa = 0.1;
% nu = 2;
% r = matern_covariance(dist, sigma2, kappa, nu);
% rho_c = 2;
% x = 1;
% r = cauchy_covariance(dist, sigma2,rho_c,x);
% r = spherical_covariance(linspace(0,0.73*rho_c,50*60), sigma2,rho_c);
% plot(dist, r)
% return

%% 2.
% First use matern_covariance to create a Sigma-covariance matrix.
% and set mu=(a constant mean).
sigma2 = 2;
kappa = 0.1;
nu = 0.1;
mu = 0;

sz = [50 60];
N = prod(sz);

% Alternative way of calculating distance matrix
% [Y, X] = meshgrid(0:sz(1) - 1, 0:sz(2) - 1);
% X = X';
% Y = Y';
% X = [X(:), Y(:)]';
% D = zeros(N, N);
% for i = 0:N - 1
%     x0 = [floor(i / sz(1)); mod(i, sz(1))];
%     x = X - x0;
%     d0 = sqrt(sum(x.^2));
%     D(:, i + 1) = d0;
% end
[u1,u2] = ndgrid(1:sz(1),1:sz(2));
D = distance_matrix([u1(:), u2(:)]);
Sigma = matern_covariance(D, sigma2, kappa, nu);
disp 'Sigma calculated';

R = chol(Sigma + eye(size(Sigma)) * 1e-5); % Calculate the Cholesky factorisation
eta = mu+R'*randn(N,1); % Simulate a sample
eta_image = reshape(eta,sz); %reshape the column to an image
%imagesc(eta_image)

%%
sigma2_epsilon = 1;
y=eta+randn(N,1)*sigma2_epsilon;
z=y-mu; 
%plot(D,z*z','.k')
%hold on 
d2 = linspace(0, max(D(:)));
figure(2)
plot(d2, matern_covariance(d2, sigma2, kappa, nu), 'r')
hold on
Kmax = 20;
Dmax = max(D(:));
[rhat,s2hat,m,n,d] = covest_nonparametric(D,z,Kmax,Dmax);
plot(d, rhat, 'b');

par_real=[sigma2; kappa; nu; sigma2_epsilon]
par=covest_ls(rhat, s2hat, m, n, d)

%%
p=0.8;
I_obs = (rand(sz)<=p);
%add nugget to the covariance matrix
Sigma_yy = Sigma + sigma_epsilon^2*eye(size(Sigma));
%and divide into observed/unobserved
Sigma_uu = Sigma_yy(~I_obs, ~I_obs);
Sigma_uo = Sigma_yy(~I_obs, I_obs);
Sigma_oo = Sigma_yy(I_obs, I_obs);
y_o = y(I_obs);
y_u = y(~I_obs);

% X = ones(prod(sz),1);
% X_u = X(~I_obs);
% X_o = X(I_obs);

y_rec = nan(sz);
y_rec(I_obs) = y_o;
y_rec(~I_obs) = mu + Sigma_uo / Sigma_oo * (y_o - mu);

V_pred = Sigma_uu - Sigma_uo / Sigma_oo * Sigma_uo';

mean(y_rec(~I_obs) - y(~I_obs))
mean(diag(V_pred))