load T_lund.mat

% plot(T_lund(:,1), T_lund(:,2))
% datetick

t = T_lund(:,1);  Y = T_lund(:,2);  n = length(Y);
X = [ones(n,1) sin(2*pi*t/365) cos(2*pi*t/365)];
beta = regress(Y, X);
eta = Y-X*beta;
% plot(t, Y, 'b', t, X*beta, 'r');
% datetick

X2 = eta(1:end-1);
Y2 = eta(2:end);
alpha = regress(Y2, X2)
e = eta(2:end) - alpha * eta(1:end-1);
sigma = var(e)
plot(t(2:end), e, 'b', t, eta, 'r');
datetick

for i = 1:25
  t2 = t + i
  X = [ones(n,1) sin(2*pi*t/365) cos(2*pi*t/365)];
  future = X * beta + alpha^i * eta
end
