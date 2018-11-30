%% Looking at GMRF
m = 100; n = 100; sz = [m, n];
kappa2 = 1e-4;
q1 = [0, -1, 0;
     -1, 4 + kappa2, -1;
     0, -1, 0];
q2 = [0, -0.1, -10;
     -0.1, 20.4 + kappa2, -0.1;
     -10, -0.1, 0];
q3 = [0, -1, -2;
     -1, 8 + kappa2, -1;
     -2, -1, 0];
Q = gmrfprec(sz, q1);
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
if false
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
end

%% Titan image (beware of runtime)
titan = double(imread('titan.jpg')) / 255;

Q = gmrfprec(size(titan), q1);
Qbeta = 1e-6 * speye(1);  %use SPEYE to create a SPARSE-identity-matrix
Qall = blkdiag(Q, Qbeta);

% miss = 0.5;
for miss = linspace(0.7, 0.95, 4)
  known=(rand(prod(size(titan)), 1)>miss);
  titan_corrupted = titan;
  titan_corrupted(~known) = 0;

  figure();
  title(miss);
  subplot(1,3,1);
  imagesc(titan);
  subplot(1,3,2);
  imagesc(titan_corrupted);

  Y = titan(known);

  A = sparse(1:length(Y), find(known), 1, length(Y), numel(titan));
  Aall = [A, ones(length(Y), 1)];

  Qeps = 1e5 * speye(length(Y));
  Qxy = Qall + Aall' * Qeps * Aall;

  p = amd(Qxy);
  Qxy = Qxy(p, p);
  Aall = Aall(:, p);

  Exy = Qxy \ Aall' * Qeps * Y;
  Exy(p) = Exy;
  Ezy = [speye(size(Q, 1)), ones(size(Q, 1), 1)] * Exy;
  subplot(1,3,3)
  imagesc(reshape(Ezy, size(titan)));
end
