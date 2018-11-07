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

