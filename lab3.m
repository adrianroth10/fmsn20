%% Looking at GMRF
m = 100; n = 100; sz = [m, n];
kappa2 = 1e-4;
q = [0, -1, 0;
     -1, 4 + kappa2, -1;
     0, -1, 0];
Q = gmrfprec(sz, q);
s = Q \ sparse(m * n / 2 + m / 2, 1, 1, size(Q, 1), 1);
% imagesc(reshape(s, m, n));
% R = chol(Q);
% x = R \ randn(size(R, 1), 1);
% figure();
% imagesc(reshape(x, m, n));

% Q2 = Q' * Q;
% figure();
% spy(Q);
% figure();
% spy(Q2)

%% Interpolation
load lab3.mat

Y = xmiss(known);

A = sparse(1:length(Y), find(known), 1, length(Y), numel(xmiss));
Aall = [A, ones(length(Y), 1)];
Qbeta = 1e-6 * speye(1);  %use SPEYE to create a SPARSE-identity-matrix
Qall = blkdiag(Q, Qbeta);

%assume a very small observation-uncertainty
Qeps = 1e5 * speye(length(Y));
%posterior precision
Qxy = Qall + Aall' * Qeps * Aall;

p = amd(Qxy);
Qxy = Qxy(p, p);
Aall = Aall(:, p);

Exy = Qxy \ Aall' * Qeps * Y;
Exy(p) = Exy;
Ezy = [speye(size(Q, 1)), ones(size(Q, 1), 1)] * Exy;
figure();
imagesc(xmiss);
figure()
imagesc(reshape(Ezy, size(xmiss)));
