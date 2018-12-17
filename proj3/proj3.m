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


n_classes = [5];
which_components = {[1,2,3]};
%which_components = {[3,4,5]};

is_beta = true;
data = beta(:,:,3:end);
proj3_pca
return;

is_beta = false;
data = img;
proj3_pca

