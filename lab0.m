x = imread('lapporten.jpg');
x = double(x)/255;
a = x(:,:,1);
b = x(:,:,2);
c = x(:,:,3);

figure(1)
imagesc(a);
%%caxis([0.5,0.9])
colormap(gray)
colorbar()
figure(2)
imagesc(b);
colorbar()
figure(3)
imagesc(c);
colorbar()

%%
x = imread('lapporten.jpg');
x = double(x)/255;
x=x/max(x(:));
img(:,:,1)=(x(:,:,1)<=0.25);
img(:,:,2)=(x(:,:,1)>0.25)&(x(:,:,1)<=0.5);
img(:,:,3)=(x(:,:,1)>0.5)&(x(:,:,1)<=0.75);
img(:,:,4)=(x(:,:,1)>0.75);
y = rgbimage(img);
imagesc(y)

%%
x = double(imread('lapporten.jpg'));
x=double(x)/255;
y=x./repmat(sum(x,3),[1,1,3]);
imagesc(y)
%%

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
for i = 1:5:25
  t2 = t + i;
  X = [ones(n,1) sin(2*pi*t2/365) cos(2*pi*t2/365)];
  future = X * beta + alpha^i * eta;
  upper = future + 1.96 * sqrt(varsum * sigma2);
  lower = future - 1.96 * sqrt(varsum * sigma2);
  figure
  hold on
  plot(t, Y, 'b', t2, X*beta, 'r');
  datetick
  plot(t2, future, 'g')
  plot(t2, upper, 'c')
  plot(t2, lower, 'c')
  xlim([t(1), t(end)])
  pause
  varsum = varsum + alpha^(2*i);
end
