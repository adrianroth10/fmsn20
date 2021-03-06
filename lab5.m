%% 1. Simple data
xc=load('lab6simple.dat');
x=reshape(xc,[65,65,2]);
% figure()
% subplot(2,2,1)
% image(rgbimage(x));
% subplot(2,2,2)
% plot(xc(:, 1), xc(:, 2), '.')
% subplot(2,2,3)
% hist2(xc, 100)
% colormap('gray')

%% 2. K-means
% [cl, theta] = kmeans(xc, 2, inf, 2);
% figure()
% imagesc(reshape(cl, 65, 65))

%% 3. GMM
% [theta, prior] = normmix_gibbs(xc, 2);
% [cl, cl_ind, p] = normmix_classify(xc, theta, prior);

%% 4. LANDSAT data
x = double(lanread('mississippi.lan')) / 255;
[y,P,Pvar] = pca(colstack(x));
% plot(Pvar)
y = y(:, 1:3);
[theta, prior] = normmix_gibbs(y, 3)
[cl, cl_ind, p] = normmix_classify(y, theta, prior);
figure()
subplot(1, 2, 1)
imagesc(reshape(y, [512, 512, 3]))
subplot(1, 2, 2)
imagesc(reshape(cl, [512, 512]))

%%
figure()
x = rgbimage(reshape(p, [512,512,3]));
subplot(1, 3, 1)
imagesc(x(:,:,1))
subplot(1, 3, 2)
imagesc(x(:,:,2))
subplot(1, 3, 3)
imagesc(x(:,:,3))

%% Parallel Gibbs sampling
N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
N3 = [0 0 0 1 1; 0 0 1 1 1; 0 1 0 1 0; 1 1 1 0 0; 1 1 0 0 0];
x = zeros(128,128,3);
iter = 100;
Plog = zeros(iter,1);
for i = 1:iter
x=mrf_sim(x,N3,0,1,1);
[~,Mz,Mf]=mrf_sim(x,N3,0,0.1,0);
Plog(i) = sum(log(Mz(logical(x))));
imagesc(x)
drawnow
end

