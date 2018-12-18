% Load data
load fmri.mat

%size of data
sz = size(img);

beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);
region_of_interest = beta(:, :, 3:end);

% imagesc(mean(region_of_interest, 3))
% colorbar;
% print('proj3/output/meanactivity.png', '-dpng');
% close

MHsigma2 = [4.7e-2, 4.5e-2, 2.8e-2, 3.5e-2, 3.3e-2     ;
            5.7e-3, 5.0e-3, 3.9e-3, 4.3e-3, 5.3e-3    ;
            4.5e-3, 4.8e-3, 3.9e-3, 3.9e-3, 3.1e-3   ;
            4.5e-3, 4.6e-3, 4.6e-3, 3.9e-3, 3.8e-3  ;
            4.3e-3, 3.8e-3, 4.2e-3, 2.8e-3, 2.8e-3];

which_components = {[1], [1,2], [1,2,3],[1,2,3,4], [1,2,3,4,5]};
n_classes = [3];

data = region_of_interest;
proj3_pca(sz, which_components, n_classes, data, MHsigma2);
