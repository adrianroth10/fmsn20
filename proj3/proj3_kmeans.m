[cl, theta] = kmeans(y_all, nc);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/kmeans_', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta), '.png'], '-dpng');
close;
