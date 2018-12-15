[theta, prior] = normmix_gibbs(y_all_stacked, nc);
[cl, cl_ind, p] = normmix_classify(y_all_stacked, theta, prior);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/gmm_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta), '.png'], '-dpng');
close;
