function proj3_gmm(y_all_stacked, sz, is_beta, str_components, nc)

[theta, prior] = normmix_gibbs(y_all_stacked, nc);
[cl, cl_ind, p] = normmix_classify(y_all_stacked, theta, prior);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/gmm_', num2str(is_beta), '_', str_components, '_', num2str(nc), '.png'], '-dpng');
close;

end
