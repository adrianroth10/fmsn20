%% Loading data and setting variables
load swissRainfall.mat

swissGrid = [swissElevation(:) swissX(:) swissY(:)];
notnanim = ~isnan(swissX);
swissGrid = swissGrid( ~isnan(swissGrid(:,1)),:);
Y = swissRain(swissRain(:,5)==0,:);
Y_valid = swissRain(swissRain(:,5)==1,:);
D = distance_matrix([Y(:, 3), Y(:, 4)]);
sz = size(swissElevation);
n = size(Y, 1);

transform = 'sqrt';
if strcmp(transform, 'sqrt')
  y = sqrt(Y(:, 1));
elseif strcmp(transform, 'log')
  y = log(Y(:, 1)+1);
end


I_valid = logical([zeros(size(Y, 1), 1); ones(size(Y_valid, 1), 1); zeros(size(swissGrid, 1), 1)]);
I_obs = logical([ones(size(Y, 1), 1); zeros(size(Y_valid, 1), 1); zeros(size(swissGrid, 1), 1)]);
coords_all = [Y(:, 3:4); Y_valid(:, 3:4); swissGrid(:, 2:3)];
elev_all = [Y(:, 2); Y_valid(:, 2); swissGrid(:, 1)];

X_u = elev_all(~I_obs, :);
X_u = [ones(size(X_u, 1), 1), X_u, X_u.^2];
X_k = elev_all(I_obs, :);
X_k = [ones(size(X_k, 1), 1), X_k, X_k.^2];

D_all = distance_matrix(coords_all);


%% Interpolation using only regression term

% Estimation of second degree polynomial for regression term
X = [ones(length(y), 1), Y(:, 2), Y(:, 2).^2];
beta = (X' * X) \ X' * y;

% Beta significance
z = y - X * beta;
s2 = var(z);
beta_var = s2 * inv((X' * X));
beta_confidence_interval = [beta - 1.96 * sqrt(diag(beta_var)), beta + 1.96 * sqrt(diag(beta_var))];
% figure()
% hold on
% lin_x = linspace(min(Y(:, 2)), max(Y(:, 2)));
% plot(Y(:, 2), y, 'b*')
% plot(lin_x, polyval(flip(beta'), lin_x), 'r')
% xlabel('Elevation (km)', 'interpreter', 'latex')
% ylabel('$\sqrt{\textnormal{precipitation}}$ ($\sqrt{\textnormal{mm}}$)', 'interpreter', 'latex')
% legend('Raw data', 'Elevation polynomial estimation')


% Interpolation
y_rec_reg = nan(size(I_obs));
y_rec_reg(I_obs) = y_k;
y_rec_reg(~I_obs) = X_u * beta;
y_rec_reg(y_rec_reg < 0) = 0;
if strcmp(transform, 'sqrt')
  y_rec_reg = y_rec_reg.^2;
elseif strcmp(transform, 'log')
  y_rec_reg = exp(y_rec_reg)-1;
end
V_reg = s2 + sum((X_u*beta_var).*X_u,2);

reg_validation_rms_error = sqrt(mean(V_reg(I_valid) .* (y_rec_reg(I_valid) - Y_valid(:, 1)).^2))

% Plotting results
rain_rec_reg = nan(sz);
rain_rec_reg(notnanim) = y_rec_reg((size(Y, 1) + size(Y_valid, 1) + 1):end);
figure()
hold on
imagesc([0 max(swissX(:))], [0 max(swissY(:))], rain_rec_reg, ...
        'alphadata', ~isnan(rain_rec_reg))
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(Y(:,3), Y(:,4), 20, Y(:,1), 'filled','markeredgecolor','r')
scatter(Y_valid(:,3), Y_valid(:,4), 20, Y_valid(:,1), 'filled','markeredgecolor','g')
xlabel('X-distance (km)');
ylabel('Y-distance (km)');
c = colorbar;
axis xy tight; hold off; c; ylabel(c, 'Rain (mm)');

figure()
hold on
plot(Y_valid(:, 1), '*')
errorbar(y_rec_reg(I_valid), 1.96 * sqrt(V_reg(I_valid)), '.')
xlim([0.5, 10.5])
xlabel('Indicies of validation data')
ylabel('Rain (mm)')
legend('Real values', 'Estimated values')




%% Is the residual covariance significant?

Dmax = max(D(:));
Kmax = 50;
[rhat, s2hat, m, n, d] = covest_nonparametric(D, z, Kmax, Dmax);
N_boot = 100;
rhats = nan(length(rhat), N_boot);
for i = 1:N_boot
  z2 = z(randperm(length(z)));
  [rhat_temp, ~, ~, ~, ~] = covest_nonparametric(D, z2, Kmax, Dmax);
  rhats(:, i) = rhat_temp;
end
rhats_sorted = sort(rhats, 2);
boot_min = rhats_sorted(:, 5);
boot_max = rhats_sorted(:, 95);
% figure();
% hold on;
% xlabel('Distance (km)')
% ylabel('Covariance')
% plot(d, boot_min, '--r')
% plot(d, boot_max, '--r')
% plot(d, rhat, 'b')
% return





%% Interpolation using both regression and covariance

% Estimating covariance matrix/field
covf = 'matern';
par_fixed = [0, 0, 0, 0];
[par, beta_refined] = covest_ml(D, y, covf, par_fixed, X);
%par = covest_ml(D, y, covf, par_fixed, X, 'reml', );

% plot(d, rhat, 'b');
% hold on
% plot(d, matern_covariance(d, par(1), par(2), par(3)));

% Beta refined significance
Sigma = matern_covariance(D, par(1), par(2), par(3));
beta_refined_var = inv((X' / Sigma * X));
beta_refined_confidence_interval = [beta_refined - 1.96 * sqrt(diag(beta_refined_var)),...
                                    beta_refined + 1.96 * sqrt(diag(beta_refined_var))];

% Interpolation
Sigma_all = matern_covariance(D_all, par(1), par(2), par(3));
sigma2_epsilon = par(4);

Sigma_yy = Sigma_all + sigma2_epsilon * eye(size(Sigma_all));
Sigma_uu = Sigma_yy(~I_obs, ~I_obs);
Sigma_uk = Sigma_yy(~I_obs, I_obs);
Sigma_kk = Sigma_yy(I_obs, I_obs);
y_k = y;

y_rec_cov = nan(size(I_obs));
y_rec_cov(I_obs) = y_k;
y_rec_cov(~I_obs) = X_u * beta_refined + Sigma_uk / Sigma_kk * (y_k - X_k * beta_refined);
y_rec_cov(y_rec_cov < 0) = 0;
if strcmp(transform, 'sqrt')
  y_rec_cov = y_rec_cov.^2;
elseif strcmp(transform, 'log')
  y_rec_cov = exp(y_rec_cov)-1;
end
V_cov_1 = diag(Sigma_uu) - sum((Sigma_uk / Sigma_kk) .* Sigma_uk, 2);
quad_coeff = (X_u' - X_k' / Sigma_kk * Sigma_uk');
V_cov_2 = sum((quad_coeff' * beta_refined_var) .* quad_coeff', 2);
V_cov = V_cov_1 + V_cov_2;

cov_validation_rms_error = sqrt(mean(V_cov(I_valid) .* (y_rec_cov(I_valid) - Y_valid(:, 1)).^2))


% Plotting results
rain_rec_cov = nan(sz);
rain_rec_cov(notnanim) = y_rec_cov((size(Y, 1) + size(Y_valid, 1) + 1):end);
figure()
hold on
imagesc([0 max(swissX(:))], [0 max(swissY(:))], rain_rec_cov, ...
        'alphadata', ~isnan(rain_rec_cov))
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(Y(:,3), Y(:,4), 20, Y(:,1), 'filled','markeredgecolor','r')
scatter(Y_valid(:,3), Y_valid(:,4), 20, Y_valid(:,1), 'filled','markeredgecolor','g')
xlabel('X-distance (km)');
ylabel('Y-distance (km)');
c = colorbar;
axis xy tight; hold off; c; ylabel(c, 'Rain (mm)');


figure()
hold on
plot(Y_valid(:, 1), '*')
errorbar(y_rec_cov(I_valid), 1.96 * sqrt(V_cov(I_valid)), '.')
xlim([0.5, 10.5])
xlabel('Indicies of validation data')
ylabel('Rain (mm)')
legend('Real values', 'Estimated values')
