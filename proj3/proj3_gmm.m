function proj3_gmm(y_all_stacked, sz, str_components, nc)

rng(0)

[theta, prior] = normmix_gibbs(y_all_stacked, nc);
[cl, cl_ind, p] = normmix_classify(y_all_stacked, theta, prior);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/gmm_', str_components, '_', num2str(nc), '.png'], '-dpng');
close;

end
