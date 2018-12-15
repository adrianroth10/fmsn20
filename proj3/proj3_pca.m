[y_data, ~, P_data] = pca(colstack(data));

for i_component = 1:length(which_components)
  components = which_components{i_component};
  for nc = n_classes
    y_all_stacked = y_data(:, :, components);
    y_all = reshape(y_all_stacked, sz(1), sz(2), []);

    proj3_kmeans;
    proj3_gmm;
    %proj3_mrf;
  end
end

