function proj3_mrf_neighbours(y_all, y_all_stacked, sz, neighbours, i_component, str_components, nc, i_neighbours, single_beta, MHsigma2)
rng(0)

iterations = 1000;
burn_in = 300;

[theta, prior] = normmix_gibbs(y_all_stacked, nc);
[~, cl_ind] = normmix_classify(y_all_stacked, theta, prior);

% change to random field
zmat = reshape(cl_ind, sz(1), sz(2), []);
zsum = zeros(sz(1), sz(2), nc);
alpha = ones(iterations + 1, nc - 1) / nc;
if single_beta
    beta = zeros(iterations + 1, 1);
    beta_prior = 1e-1;
else
    beta = zeros(iterations + 1, nc);
    beta_prior = 1e-1 * eye(nc);
end
Plog = zeros(iterations, 1);
acc = 0;

for iter = 1:iterations
    alpha_post = mrf_gaussian_post([0, alpha(iter, :)], theta, y_all);
    zmat = mrf_sim(zmat, neighbours, alpha_post, beta(iter, :), 1);
    [~, Mz] = mrf_sim(zmat, neighbours, alpha_post, beta(iter, :), 0);
    Plog(iter) = sum(log(Mz(logical(zmat))));
    for k = 1:nc
        if sum(sum(zmat(:,:,k)))>=2
            ind = find(zmat(:, :, k));
            [mu, Sigma] = gibbs_mu_sigma(y_all_stacked(ind, :));
            theta{k}.mu = mu;
            theta{k}.Sigma = Sigma;
        end
    end
    [alpha(iter + 1, :), beta(iter + 1, :), acc_tmp] = ...
                  gibbs_alpha_beta(alpha(iter, :), beta(iter, :), zmat, neighbours, beta_prior, MHsigma2(i_component, nc - 1, i_neighbours) / (2 - single_beta));
    acc = acc + acc_tmp;
    if iter > burn_in
        zsum = zsum + zmat;
    end
end
acc = acc / iterations;
fprintf(1, 'acceptance rate = %.4f\n', acc);

matrix2latex(['proj3/output/acceptance_rate_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '_', num2str(single_beta), '.tex'], acc);
[zest, xest] = max(zsum / iterations, [], 3);

figure();
imagesc(xest);
print(['proj3/output/mrf_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '_', num2str(single_beta), '.png'], '-dpng');
close

figure();
plot(Plog);
xlabel('Gibbs iteration')
ylabel('pseudo log likelihood of field')
print(['proj3/output/mrf_plog_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '_', num2str(single_beta), '.png'], '-dpng');
close

figure();
plot(alpha(1:iterations, :));
xlabel('Gibbs iteration')
ylabel('\alpha values')
print(['proj3/output/mrf_alpha_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '_', num2str(single_beta), '.png'], '-dpng');
close

figure();
plot(beta(1:iterations, :));
xlabel('Gibbs iteration')
ylabel('\beta values')
print(['proj3/output/mrf_beta_', str_components, '_', num2str(nc), '_', num2str(i_neighbours), '_', num2str(single_beta), '.png'], '-dpng');
close

end
