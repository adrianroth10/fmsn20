function proj3_kmeans(y_all_stacked, sz, str_components, nc)

rng(0);

[cl, theta] = kmeans(y_all_stacked, nc);

figure();
imagesc(reshape(cl, sz(1:2)));
print(['proj3/output/kmeans_', str_components, '_', num2str(nc), '.png'], '-dpng');
close;

end
