% Load data
load fmri.mat

%size of data
sz = size(img);

beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);


n_classes = [2, 3, 4, 5];
is_beta = true;
data = beta;
which_components = {[1], [1, 2], [1,2,3], [1,4,5]};
proj3_pca
return;

is_beta = false;
data = img;
which_components = {[1], [1, 2], [1,2,3], [1,4,5]};
proj3_pca
