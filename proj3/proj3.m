% Set seed
rng(0)

% Load data
load fmri.mat

%size of data
sz = size(img);

beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);

n_classes = [5];
is_beta = true;
data = beta;
which_components = {[1]};
proj3_pca


is_beta = false;
data = img;
which_components = { [2], [2,3], [2,3 4]};
proj3_pca
