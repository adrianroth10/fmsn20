load swissRainfall.mat

swissGrid = [swissElevation(:) swissX(:) swissY(:)];
notnanim = ~isnan(swissX);
swissGrid = swissGrid( ~isnan(swissGrid(:,1)),:);
Y = swissRain(swissRain(:,5)==0,:);
Yvalid = swissRain(swissRain(:,5)==1,:);

sz = size(swissElevation);
n = size(Y, 1);
X = [ones(n, 1), Y(:, 2), Y(:, 2).^2];
y = Y(:, 1);

% Estimation of second degree polynomial for regression term
beta = (X' * X) \ X' * y;
% plot(Y(:, 2), Y(:, 1) - X * beta, '*')
z = y - X * beta;

% Estimating covariance matrix/field
D = distance_matrix([Y(:, 3), Y(:, 4)]);

covf = 'matern';
par_fixed = [0, 0, 0.9, 0];
Kmax = 50;
Dmax = max(D(:));
[rhat, s2hat, m, n, d] = covest_nonparametric(D, z, Kmax, Dmax);
par = covest_ls(rhat, s2hat, m, n, d, covf, par_fixed);

% Refinement step
Sigma = matern_covariance(D, par(1), par(2), par(3));
beta_refined = (X' / Sigma * X) \ X' / Sigma * y;
z_refined = y - X * beta_refined;
[rhat, s2hat, m, n, d] = covest_nonparametric(D, z_refined, Kmax, Dmax);
par = covest_ls(rhat, s2hat, m, n, d, covf, par_fixed);

% plot(d, rhat, 'b');
% hold on
% plot(d, matern_covariance(d, par(1), par(2), par(3)));

% Interpolation
I_obs = logical([ones(size(Y, 1), 1); zeros(size(Yvalid, 1), 1); zeros(size(swissGrid, 1), 1)]);
coords_all = [Y(:, 3:4); Yvalid(:, 3:4); swissGrid(:, 2:3)];
elev_all = [Y(:, 2); Yvalid(:, 2); swissGrid(:, 1)];


D_all = distance_matrix(coords_all);
Sigma_all = matern_covariance(D_all, par(1), par(2), par(3));

% add nugget to the covariance matrix
Sigma_yy = Sigma_all + par(4)*eye(size(Sigma_all));
% and divide into observed/unobserved
Sigma_uu = Sigma_yy(~I_obs, ~I_obs);
Sigma_uk = Sigma_yy(~I_obs, I_obs);
Sigma_kk = Sigma_yy(I_obs, I_obs);
y_k = y;

X_u = elev_all(~I_obs, :);
X_u = [ones(size(X_u, 1), 1), X_u, X_u.^2];
X_k = elev_all(I_obs, :);
X_k = [ones(size(X_k, 1), 1), X_k, X_k.^2];

y_rec = nan(size(I_obs));
y_rec(I_obs) = Y(:, 1);
y_rec(~I_obs) = X_u * beta_refined + Sigma_uk / Sigma_kk * (y_k - X_k * beta_refined);

rain_rec = nan(sz);
rain_rec(notnanim) = y_rec((size(Y, 1) + size(Yvalid, 1) + 1):end);
imagesc(rain_rec)
axis xy tight; hold off; colorbar
