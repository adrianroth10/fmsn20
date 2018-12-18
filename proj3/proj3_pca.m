function proj3_pca(sz, which_components, n_classes, is_beta, data, MHsigma2)
[y_data, ~, P_data] = pca(colstack(data));
% plot(P_data)
% print(['proj3/output/pdata_', num2str(is_beta),'.png'], '-dpng');
% close;
for i_component = 1:length(which_components)
  components = which_components{i_component};
  str_components = strrep(num2str(components), ' ', '');
  for nc = n_classes
    y_all_stacked = y_data(:, components);
    y_all = reshape(y_all_stacked, sz(1), sz(2), []);

    %proj3_kmeans(y_all_stacked, sz, is_beta, str_components, nc);
    %proj3_gmm(y_all_stacked, sz, is_beta, str_components, nc);
    proj3_mrf(y_all, y_all_stacked, sz, is_beta, length(components), str_components, nc, MHsigma2) ;

  end
end

end
