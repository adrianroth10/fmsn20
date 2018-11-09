%% 1.

%dist = 0:80;
%sigma2 = 1;
%kappa = 0.1;
%nu = 2;
%r = matern_covariance(dist, sigma2, kappa, nu);
%rho_c = 2;
%x = 1;
%r = cauchy_covariance(dist, sigma2,rho_c,x);
%r = spherical_covariance(linspace(0,0.73*rho_c,50*60), sigma2,rho_c);
%plot(dist, r)
%return

%% 2.
% First use matern_covariance to create a Sigma-covariance matrix.
% and set mu=(a constant mean).
sigma2 = 1;
kappa = 0.1;
nu = 1;
mu = 0;

sz = [50 52];
N = prod(sz);
[Y, X] = meshgrid(0:sz(1) - 1, 0:sz(2) - 1);
X = [X(:), Y(:)]';

d = zeros(N, N);
for i = 0:N - 1
    x0 = [mod(i, sz(2)); floor(i / sz(2))];
    x = X - x0;
    d0 = sqrt(sum(x.^2));
    d(:, i + 1) = d0;
end
Sigma = matern_covariance(d, sigma2, kappa, nu);
disp 'Sigma calculated';

R = chol(Sigma + eye(size(Sigma)) * 1e-5); % Calculate the Cholesky factorisation
eta = mu+R'*randn(N,1); % Simulate a sample
eta_image = reshape(eta,sz); %reshape the column to an image
imagesc(eta_image)
