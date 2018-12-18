% Set seed
rng(0)

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

MHsigma2 = [0, 4.5e-2, 0, 0, 0;
            0, 4.5e-3, 0, 0, 0;
            0, 3.5e-3, 0, 0, 0;
            0, 2.5e-3, 0, 0, 0;
            0, 2.0e-3, 0, 0, 0;
            0, 1.5e-3, 0, 0, 0];

which_components = {[1], [1,2], [1,2,3], [1,2,3,4], [1,2,3,4,5]};
% which_components = {[1], [1,2], [1,2,3]};
n_classes = [3];

is_beta = true;
data = region_of_interest;
proj3_pca(sz, which_components, n_classes, is_beta, data, MHsigma2);
return;

is_beta = false;
data = img;
proj3_pca

