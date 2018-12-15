[y_data, ~, P_data] = pca(colstack(data));

for i_component = 1:length(which_components)
  component = which_components{i_component};
  for nc = n_classes
    y_all = y_data(:, component);

    proj3_kmeans;
    proj3_gmm;
    %proj3_mrf;
  end
end

