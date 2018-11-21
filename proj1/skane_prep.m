%% Finding border
im_border = imread('skane.png') == 0;
im_border_orig = im_border;
s = size(im_border);
[r_all, c_all] = find(im_border);
start = [r_all(1), c_all(1)];
curr = start;
n = 1;
coords = zeros(sum(im_border(:)), 2);
coords(1, :) = curr;
while 1
  r = curr(1); c = curr(2);
  rgz = r - 1 > 0;
  rlm = r + 1 < s(1);
  cgz = c - 1 > 0;
  clm = c + 1 < s(2);
  if rgz && im_border(r - 1, c)
    curr = [r - 1, c];
  elseif cgz && im_border(r, c - 1)
    curr = [r, c - 1];
  elseif rlm && im_border(r + 1, c)
    curr = [r + 1, c];
  elseif clm && im_border(r, c + 1)
    curr = [r, c + 1];
  elseif rgz && cgz && im_border(r - 1, c - 1)
    curr = [r - 1, c - 1];
  elseif rgz && clm && im_border(r - 1, c + 1)
    curr = [r - 1, c + 1];
  elseif rlm && cgz && im_border(r + 1, c - 1)
    curr = [r + 1, c - 1];
  elseif rlm && rlm && im_border(r + 1, c + 1)
    curr = [r + 1, c + 1];
  else
    break;
  end
  im_border(curr(1), curr(2)) = 0;
  n = n + 1;
  coords(n, :) = curr;
end
coords(n + 1, :) = start;
coords = coords(1:n + 1, :);

%% Estimating scale and coordinates
width_scale = (14.57 - 12.46) / (max(r_all) - min(r_all));
width_const = 14.57 - width_scale * max(r_all);
height_scale = (56.52 - 55.34) / (max(c_all) - min(c_all));
height_const = 56.52 - height_scale * max(c_all);

skaneBorder = fliplr(coords);
skaneBorder(:, 1) = skaneBorder(:, 1) * width_scale + width_const;
skaneBorder(:, 2) = skaneBorder(:, 2) * height_scale + height_const;

notnanim = imfill(im_border_orig, 'holes');
% notnanim =

[X, Y] = meshgrid(1:s(1), 1:s(2));
X = X * width_scale * 1.06 + width_const - 0.04;
Y = Y * height_scale * 1.06 + height_const - 0.03;
% small time tweaking for perfection

skaneElevation = nan(s);
skaneElevation(notnanim) = 0; % skåne is flat
skaneX = nan(s);
skaneX(notnanim) = X(notnanim);
skaneY = nan(s);
skaneY(notnanim) = Y(notnanim);

%% Measured values with validation data
data = csvread('proj1/skane_2018-08-12_data.csv', 0, 1);
n_points = size(data, 1);

skaneRain = zeros(n_points, 5);
skaneRain(:, 1) = data(:, 4);
skaneRain(:, 2) = zeros(n_points, 1);
skaneRain(:, 3) = data(:, 3);
skaneRain(:, 4) = data(:, 2);
skaneRain(:, 5) = randperm(n_points) <= 8;

save('proj1/skaneRainfall.mat', 'skaneRain', 'skaneElevation', 'skaneX', 'skaneY', 'skaneBorder');

%% Testplot
I_valid = skaneRain(:, 5) == 1;
figure(); hold on;
imagesc([min(skaneX(:)) max(skaneX(:))], [min(skaneY(:)) max(skaneY(:))], skaneElevation, ...
        'alphadata', ~isnan(skaneElevation))
plot(skaneBorder(:,1), skaneBorder(:,2),'k')
scatter(skaneRain(~I_valid, 3), skaneRain(~I_valid, 4), 20, skaneRain(~I_valid, 1), 'filled','markeredgecolor','r')
scatter(skaneRain(I_valid, 3), skaneRain(I_valid, 4), 20, skaneRain(I_valid, 1), 'filled','markeredgecolor','g')
