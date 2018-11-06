function marks=bioid_landmarks(indices,thepath)
% BIOID_LANDMARKS Read face landmark data from the BioID-dataset.
%
% marks=bioid_landmarks(indices)
% marks=bioid_landmarks(indices,thepath)
%
% indices: vector of image indices
% marks  : p-by-2-by-n matrix of landmarks, p landmarks of n faces
% thepath: The path to the *.pts-files.  By default, we assume
%          that the files are located in 
%          ./data/protected/face/bioid/points_20/,
%          which is the standard on the EFD computer network
%          If the fucntion can find an image file it will revert to the
%          default path for a second attempt before failing.
%
% Data:
%   http://www.humanscan.de/support/downloads/facedb.php
%
% See also: bioid_images

% Copyright (c) 2004-2005 Finn Lindgren
% $Id: bioid_landmarks.m 3006 2006-10-05 14:25:17Z johanl $

if (nargin<2), thepath = []; end
if (isempty(thepath)),
  layout = 1;
  [p,n,e]=fileparts(which(mfilename));
  thepath = fullfile(p,filesep,'data',filesep,'protected',filesep, ...
		     'face',filesep,'bioid',filesep, ...
		     'points_20',filesep); 
else
  layout = 0;
end

if (length(indices)<1)
  marks = zeros(20,2,0);
  return;
end

marks = zeros(20,2,length(indices));
for k=1:length(indices)
  index = indices(k);
  filename = sprintf('%sbioid_%04i.pts',thepath,index);
  fid = fopen(filename,'r');
  if (fid<0)
    % Try other directory layout.
    if (layout==0)
      [p,n,e]=fileparts(which(mfilename));
      thepath = fullfile(p,filesep,'data',filesep,'protected',...
                         filesep,'face',filesep,'bioid',...
                         filesep,'points_20',filesep);
      filename = sprintf('%sbioid_%04i.pts',thepath,index);
      fid = fopen(filename,'r');
    end
    if (fid<0) % If not found anywhere
      warning(sprintf('Could not open file: %s',filename));
      continue;
    end
  end
  fgetl(fid); % version: 1
  fgetl(fid); % n_points: 20
  fgetl(fid); % {
  [result,count] = fscanf(fid,'%f'); % 159.128 108.541, on each row
  % }
  if (count ~= 40)
    warning(sprintf('Trouble parsing: %s',filename));
    continue;
  end
  fclose(fid);
  marks(:,:,k) = reshape(result,[2,20])';
end
