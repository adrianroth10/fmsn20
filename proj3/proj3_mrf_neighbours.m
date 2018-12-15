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
alpha = zeros(nc, iterations + 1);
beta = zeros(iterations + 1);
Plog = zeros(iterations);
acc = 0;

for iter = 1:iterations
    alpha_post = mrf_gaussian_post(alpha(:, iter)', theta, y_all);
    zmat = mrf_sim(zmat, neighbours, alpha_post, beta(iter), 1);
    [~, Mz] = mrf_sim(zmat, neighbours, alpha_post, beta(iter), 0);
    Plog(iter) = sum(log(Mz(logical(zmat))));
    for k = 1:nc
        [mu, Sigma] = gibbs_mu_sigma(y_all(logical(zmat(:, :, k))));
        theta{k}.mu = mu;
        theta{k}.Sigma = Sigma;
    end
    [alpha(2:end, iter + 1), beta(iter + 1), acc_tmp] = ...
                  gibbs_alpha_beta(alpha(2:end, iter)', beta(iter), zmat, neighbours, 10, 1e-2);
    acc = acc + acc_tmp;
    if iter > burn_in
        zsum = zsum + zmat;
    end
end
acc = acc / iterations
[zest, xest] = max(zsum / iterations, [], 3);

figure();
imagesc(xest);
print(['proj3/output/mrf_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(Plog);
print(['proj3/output/mrf_plog_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(alpha);
print(['proj3/output/mrf_alpha_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbours), '.png'], '-dpng');
close;

figure();
plot(beta);
print(['proj3/output/mrf_beta_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbours), '.png'], '-dpng');
close;
