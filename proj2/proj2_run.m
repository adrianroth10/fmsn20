% type: 1 (CAR), 2 (SAR), 3 (OSC SAR)
if type == 1 || type == 2
  par_init = [0 0];
else
  par_init = [0 0 0];
end

global x_mode;
x_mode = [];
par = fminsearch( @(x) GMRF_negloglike(x, Y(I), Atilde(I,:), C, ...
  G, G2, qbeta, type), par_init);

E_xy = Agrid * x_mode(1:end - Nbeta);
E_beta = x_mode(end - Nbeta + 1:end);
E_Bbeta = Bgrid * E_beta;
E_zy = E_xy + E_Bbeta;
E_out = exp(E_zy);

%reuse taylor expansion to compute posterior precision
tau = exp(par(1));
kappa2 = exp(par(2));
if type == 1
  Q_x = tau*(kappa2*C + G);
elseif type == 2
  Q_x = tau*(kappa2^2*C + 2*kappa2*G + G2);
else
  gamma = (exp(par(3)) - 1) / (exp(par(3)) + 1);
  Q_x = tau*(kappa2^2*C + 2*gamma*kappa2*G + G2);
end

Qtilde = blkdiag(Q_x, qbeta*speye(Nbeta));

[~, ~, Q_xy] = GMRF_taylor(x_mode, Y(I), Atilde(I, :), Qtilde);

%% Variance of beta
e = [zeros(size(Q_xy,1)-Nbeta, Nbeta); eye(Nbeta)];
V_beta0 = e'*(Q_xy\e);
E_beta;
beta_int_size = sqrt(diag(V_beta0)) * 1.96;
for i = 1:Nbeta
  fprintf(1, 'Beta %d: %11.4f +- %.4f\n', i, E_beta(i), beta_int_size(i));
end

%1000 samples from the approximate posterior
Rxy = chol(Q_xy);
x_samp = bsxfun(@plus, x_mode, Rxy \ randn(size(Rxy,1), 1000));
Vx  = 1 / (size(x_samp,2) - Nbeta - length(par_init)) * sum(x_samp(1:end-Nbeta,:).^2,2);
Vbeta = sum((Bgrid * V_beta0) .* Bgrid, 2);
Vzy = Vx + Vbeta;
Vout = exp(Vzy);

rms_error = sqrt(mean((E_out(Ivalid) - Y(Ivalid)).^2));
weighted_rms_error = sqrt(mean(Vout(Ivalid) .* (E_out(Ivalid) - Y(Ivalid)).^2));
fval = GMRF_negloglike(par, Y(I), Atilde(I, :), C, G, G2, qbeta, type);
fprintf(1, 'rms_error: %.4f\nweighted_rms_error: %.4e\nfval: %.4e\n', rms_error, weighted_rms_error, fval);

%% Plotting
% Complete
figure()
imagesc(reshape(E_out, sz))
title('Estimated counts')
colorbar
return
% Mean
figure('Position', [100, 100, 1000, 400]);
sgtitle('Mean component')
subplot(1,2,1)
imagesc(reshape(E_Bbeta, sz))
title('Estimated counts')
colorbar
subplot(1,2,2)
imagesc(reshape(sqrt(Vbeta),sz))
title('Standard deviation')
colorbar
% Spatial
figure('Position', [100, 100, 1000, 400]);
sgtitle('Spatial smooth component')
subplot(1,2,1)
title('Spatial component')
imagesc(reshape(E_xy, sz))
title('Estimated counts')
colorbar
subplot(1,2,2)
imagesc(reshape(sqrt(Vx),sz))
title('Standard deviation')
colorbar
% log complete
figure('Position', [100, 100, 1000, 400]);
sgtitle('Complete log intensity')
title('Log intensity of complete component')
subplot(1,2,1)
imagesc(reshape(E_zy, sz))
title('Estimated counts')
colorbar
subplot(1,2,2)
imagesc(reshape(sqrt(Vzy),sz))
title('Standard deviation')
colorbar
% Plot variance to distance
figure()
plot(E_zy, Vzy, '*')
