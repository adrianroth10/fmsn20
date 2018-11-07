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
sigma2 = var(e)
% plot(t(2:end), e, 'b', t, eta, 'r');
% datetick

varsum = 1;
for i = 1:5
  t2 = t + i;
  X = [ones(n,1) sin(2*pi*t/365) cos(2*pi*t/365)];
  future = X * beta + alpha^i * eta;
  upper = future + 1.96 * sqrt(varsum * sigma2);
  lower = future - 1.96 * sqrt(varsum * sigma2);
  figure
  hold on
  plot(t, Y, 'b', t, X*beta, 'r');
  xlim([t(1), t(end)])
  datetick
  plot(t2, future, 'g')
  plot(t2, upper, 'c')
  plot(t2, lower, 'c')
  pause
  varsum = varsum + alpha^(2*i);
end
