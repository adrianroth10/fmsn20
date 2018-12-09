%load data
load HA2_forest
%size of the grid
sz = size(bei_counts);
%observations
Y = bei_counts(:);
% data
rng(0)  % setting seed for predictable random sequence
I = ~isnan(Y);
ind_i = find(I);
Ivalid = logical(zeros(size(I)));
Ivalid(ind_i(rand(size(ind_i)) < 0.1)) = true;
I(Ivalid) = false;

%create Q-matrix
[u1, u2] = ndgrid(1:sz(1),1:sz(2));
[C,G,G2] = matern_prec_matrices([u1(:) u2(:)]);
%mean value-vector (might not need all)
% Bgrid = [ones(prod(sz),1) bei_elev(:)];
% Bgrid = [ones(prod(sz),1) bei_grad(:)];
% % Bgrid = [ones(prod(sz),1) bei_elev(:) bei_grad(:)];
% Bgrid = [ones(prod(sz),1) bei_elev(:) bei_elev(:).^2];
% Bgrid = [ones(prod(sz),1) bei_grad(:) bei_grad(:).^2];
Bgrid = [ones(prod(sz),1) bei_elev(:) bei_elev(:).^2 bei_grad(:) bei_grad(:).^2];
Nbeta = size(Bgrid,2);
%Bgrid = [zeros(prod(sz),1)];
qbeta = 1e-6;
%and observation matrix for the grid
%Agrid = speye(prod(sz));
Agrid = speye(prod(sz));
%G2 is the most dense of the matrices, lets reorder
p = amd(G2);
%reorder precision matrices
C = C(p,p);
G = G(p,p);
G2 = G2(p,p);
%and observation matrix
Agrid = Agrid(:,p);
%create A tilde matrix
Atilde = [Agrid Bgrid];

type = 1; % car model
proj2_run;

%%
% figure()
% subplot(2,2,1)
% title('Counted data')
% imagesc(reshape(bei_counts, sz))
% colorbar
% subplot(2,2,2)
% title('Counted data')
% imagesc(reshape(Ivalid, sz))
% colorbar
% subplot(2,2,3)
% title('Elevation')
% imagesc(reshape(bei_elev,sz))
% colorbar
% subplot(2,2,4)
% title('Elevation gradient')
% imagesc(reshape(bei_grad,sz))
% colorbar
