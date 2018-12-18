function proj3_kmeans(y_all_stacked, sz, is_beta, str_components, nc)

[cl, theta] = kmeans(y_all_stacked, nc);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/kmeans_', num2str(is_beta), '_', str_components, '_', num2str(nc), '.png'], '-dpng');
close;

end
