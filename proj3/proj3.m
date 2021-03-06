% Load data
load fmri.mat

%size of data
sz = size(img);

beta = X\colstack(img)';
beta = reshape(beta', sz(1), sz(2), []);
region_of_interest = beta(:, :, 3:end);

visfile ='proj3/output/visualize1.png';
if exist(visfile, 'file') == 0
  imagesc(img(:, :, 1))
  colorbar;
  print(visfile, '-dpng');
  imagesc(img(:, :, 10))
  colorbar;
  print('proj3/output/visualize2.png', '-dpng');
  close
end

meanfile ='proj3/output/meanactivity.png';
if exist(meanfile, 'file') == 0
  imagesc(mean(region_of_interest, 3))
  colorbar;
  print(meanfile, '-dpng');
  close
end


file ='proj3/output/x.png';
if exist(file, 'file') == 0
  plot(X)
  print(file, '-dpng');
  close
end

% MHsigma2 = [4.7e-2, 4.5e-2, 2.8e-2, 3.5e-2, 3.3e-2     ;
%             5.7e-3, 5.0e-3, 3.9e-3, 4.3e-3, 5.3e-3    ;
%             4.5e-3, 4.8e-3, 3.9e-3, 3.9e-3, 3.1e-3   ;
%             4.5e-3, 4.6e-3, 4.6e-3, 3.9e-3, 3.8e-3  ;
%             4.3e-3, 3.8e-3, 4.2e-3, 2.8e-3, 2.8e-3];
% matrix2latex('proj3/output/mhsigma.tex', MHsigma2);
MHsigma2All = zeros(4, 3, 2);
MHsigma2N1 = [4.7e-2, 1.0e-2, 2.8e-2    ;
              5.7e-3, 5.0e-3, 4.1e-3   ;
              4.5e-3, 4.8e-3, 3.9e-3  ;
              4.5e-3, 4.7e-3, 4.6e-3];
matrix2latex('proj3/output/mhsigman1.tex', MHsigma2N1);

MHsigma2N2 = [1.1e-2, 0.4e-2, 0.3e-2    ;
              1.1e-3, 1.0e-3, 1.0e-3   ;
              1.0e-3, 1.0e-3, 1.0e-3  ;
              0.9e-3, 1.0e-3, 0.9e-3];
matrix2latex('proj3/output/mhsigman2.tex', MHsigma2N2);
MHsigma2All(:, :, 1) = MHsigma2N1;
MHsigma2All(:, :, 2) = MHsigma2N2;

% which_components = {[1,2,3,4]};
% n_classes = [3];
which_components = {[1], [1,2], [1,2,3], [1,2,3,4]};
n_classes = [2,3,4];

data = region_of_interest;
proj3_pca(sz, which_components, n_classes, data, MHsigma2All);
