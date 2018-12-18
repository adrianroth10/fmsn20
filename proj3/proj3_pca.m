function proj3_pca(sz, which_components, n_classes, data, MHsigma2)

[y_data, ~, P_data] = pca(colstack(data));

pdatafile ='proj3/output/pdata.png';
if exist(pdatafile, 'file') == 0
  plot(P_data)
  xlabel('Sorted eigen value #')
  ylabel('Eigen value size')
  print(pdatafile, '-dpng');
  close
end

for i_component = 1:length(which_components)
  components = which_components{i_component};
  str_components = strrep(num2str(components), ' ', '');
  for nc = n_classes
    y_all_stacked = y_data(:, components);
    y_all = reshape(y_all_stacked, sz(1), sz(2), []);

    proj3_kmeans(y_all_stacked, sz, str_components, nc);
    %proj3_gmm(y_all_stacked, sz, str_components, nc);
    %proj3_mrf(y_all, y_all_stacked, sz, length(components), str_components, nc, MHsigma2) ;

  end
end

end
