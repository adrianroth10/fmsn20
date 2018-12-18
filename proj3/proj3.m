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

MHsigma2 = [4.5e-2, 4.5e-2, 2.3e-2, 2.31e-2, 2.32e-2, 2.30e-2     ;
            4.9e-3, 4.5e-3, 2.4e-3, 2.37e-3, 2.32e-3, 2.37e-3    ;
            4.5e-3, 3.5e-3, 2.4e-3, 2.30e-3, 2.30e-3, 2.31e-3   ;
            4.5e-3, 2.5e-3, 2.3e-3, 2.34e-3, 2.30e-3, 2.31e-3  ;
            4.5e-3, 2.0e-3, 2.4e-3, 2.31e-3, 2.30e-3, 2.31e-3];
matrix2latex('proj3/output/mhsigma2.tex', MHsigma2);

which_components = {[1], [1,2], [1,2,3],[1,2,3,4], [1,2,3,4,5]};
n_classes = [2, 3, 4, 5, 6, 7];

data = region_of_interest;
proj3_pca(sz, which_components, n_classes, data, MHsigma2);
