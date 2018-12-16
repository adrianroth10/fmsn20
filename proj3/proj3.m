% Set seed
rng(0)

% Load data
load fmri.mat

%size of data
sz = size(img);

beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);


n_classes = [4 5 6 7 ];
which_components = { [1,2] [2,3] [2,3,4]};

is_beta = true;
data = beta;
proj3_pca

is_beta = false;
data = img;
proj3_pca

