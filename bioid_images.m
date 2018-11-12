function img = bioid_images(indices,thepath)
% BIOID_IMAGES Read face image data from the BioID-dataset.
%
% images=bioid_images(indices)
% images=bioid_images(indices,thepath)
%
% indices: vector of image indices
% images : a-by-b-n matrix of a-by-b images, n images
% thepath: The path to the *.mat-files.  By default, we assume
%          that the files are located in the same directory as when
%          installed on the EFD computer network, i.e.
%          ./data/protected/face/bioid/
%          If the fucntion can't find an image file it will revert to the
%          default path for a second attempt before failing.
%
% Data:
%   http://www.humanscan.de/support/downloads/facedb.php
%
% See also: bioid_landmarks

% Copyright (c) 2004 Finn Lindgren
% $Id: bioid_images.m 3006 2006-10-05 14:25:17Z johanl $

if (nargin<2), thepath = []; end
if (isempty(thepath))
  [p,n,e]=fileparts(which(mfilename));
  thepath = fullfile(p,filesep,'data',filesep,'protected',filesep, ...
		     'face',filesep,'bioid',filesep); 
end

if (length(indices)<1)
  img = zeros(286,384,0,'uint8');
  return;
end

img = zeros(286,384,length(indices),'uint8');
%figure out which mat file each image is stored in
file_index = floor( (indices-1)/100 )+1;

for k=min(file_index(:)):max(file_index(:))
  I = file_index==k;
  if( ~any(I) )
    continue;
  end;
  
  filename = sprintf('%sbioid_%04i_%04i.mat',thepath,100*(k-1)+1,100*k);
  try
    load(filename);
  catch %failed, try default path.
    [p,n,e]=fileparts(which(mfilename));
    thepath = fullfile(p,filesep,'data',filesep,'protected',filesep, ...
		       'face',filesep,'bioid',filesep); 
    filename = sprintf('%sbioid_%04i_%04i.mat',thepath,100*(k-1)+1,100*k);
    load(filename);
  end;

  ind = false(100,1);
  ind(indices(I)-100*(k-1)) = true;
  img(:,:,I) = images(:,:,ind);
end


% $$$ for k=1:length(indices)
% $$$   index = indices(k);
% $$$   i1 = index(
% $$$   filename = sprintf('%sbioid_%04i.pgm',thepath,index);
% $$$   try
% $$$     images(:,:,k) = imread(filename);
% $$$   catch %failed, try default path.
% $$$     thepath = fullfile('.',filesep,'data',filesep,'protected',filesep, ...
% $$$ 		       'face',filesep,'bioid',filesep); 
% $$$     filename = sprintf('%sbioid_%04i.pgm',thepath,index);
% $$$     images(:,:,k) = imread(filename);
% $$$   end;
% $$$ end
