iterations = 1000;
burn_in = 200;

%[theta, ~] = normmix_gibbs(y_all, nc);
theta = cell(nc,1);
for k = 1:nc
    theta{k}.mu = zeros(length(components),1);
    theta{k}.Sigma = eye(length(components));
end

zmat = zeros(sz(1), sz(2), nc);
zsum = zeros(sz(1), sz(2), nc);
alpha = zeros(iterations + 1, nc);
beta = zeros(iterations + 1);
Plog = zeros(iterations);
acc = 0;

for iter = 1:iterations
    alpha_post = mrf_gaussian_post(alpha(iter, :)', theta, y_all);
    zmat = mrf_sim(zmat, neighbours, alpha_post, beta(iter), 1);
    [~, Mz] = mrf_sim(zmat, neighbours, alpha_post, beta(iter), 0);
    Plog(iter) = sum(log(Mz(logical(zmat))));
    for k = 1:nc
        [mu, Sigma] = gibbs_mu_sigma(y_all(logical(zmat(:, :, k))));
        theta{k}.mu = mu;
        theta{k}.Sigma = Sigma;
    end
    [alpha(iter + 1, 2:end), beta(iter + 1), acc_tmp] = ...
                  gibbs_alpha_beta(alpha(iter, 2:end)', beta(iter), zmat, neighbours, 1e-1, 1e-2);
    acc = acc + acc_tmp;
    if iter > burn_in
        zsum = zsum + zmat;
    end
end
acc = acc / iterations;
fprintf(1, 'acceptance rate = %.4f\n', acc);

[zest, xest] = max(zsum / iterations, [], 3);

figure();
imagesc(xest);
print(['proj3/output/mrf_', num2str(is_beta), '_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(Plog);
xlabel('Gibbs iteration')
ylabel('pseudo log likelihood of field')
print(['proj3/output/mrf_plog_', num2str(is_beta), '_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(alpha(1:iterations, :));
xlabel('Gibbs iteration')
ylabel('\alpha values')
print(['proj3/output/mrf_alpha_', num2str(is_beta), '_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(beta(1:iterations));
xlabel('Gibbs iteration')
ylabel('\beta values')
print(['proj3/output/mrf_beta_', num2str(is_beta), '_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '.png'], '-dpng');
close;
