dist = 0:80;
sigma2 = 1;
kappa = 0.1;
nu = 4;
rho_c = 2;
x = 1;
%r = matern_covariance(dist, sigma2, kappa, nu);
%r = cauchy_covariance(dist, sigma2,rho_c,x);
%r = spherical_covariance(linspace(0,0.73*rho_c,50*60), sigma2,rho_c);
%plot( r)


% First use matern_covariance to create a Sigma-covariance matrix.
% and set mu=(a constant mean).
mu = 5;
sz = [10 20];
N = prod(sz);
r = matern_covariance(linspace(0, 60, prod(sz)), sigma2, kappa, nu);
Sigma = diag(r(1) * ones(N, 1));
for i = 2:N
    du = diag(r(i) * ones(N - i + 1, 1), i-1);
    dl = du';
    Sigma = Sigma + du + dl;
end
disp 'done';
R = chol(Sigma); % Calculate the Cholesky factorisation
eta = mu+R'*randn(60,1); % Simulate a sample
eta_image = reshape(eta,[60,60]); %reshape the column to an image
imagesc(eta_image)