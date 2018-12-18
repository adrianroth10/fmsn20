function proj3_mrf(y_all, y_all_stacked, sz, i_component, str_components, nc, MHsigma2)

% Neighbourhood
neighbours1 = [0 1 0;
               1 0 1;
               0 1 0];
neighbours2 = [1 1 1;
               1 0 1;
               1 1 1];
neighbours3 = [0 0 1 0 0
               0 1 1 1 0;
               1 1 0 1 1;
               0 1 1 1 0;
               0 0 1 0 0];
%neighbours_set = {neighbours1, neighbours2, neighbours3};
neighbours_set = {neighbours1};

for i_neighbours = 1:length(neighbours_set)
    neighbours = neighbours_set{i_neighbours};
    proj3_mrf_neighbours(y_all, y_all_stacked, sz, neighbours, i_component, str_components, nc, i_neighbours, MHsigma2);
end

end
