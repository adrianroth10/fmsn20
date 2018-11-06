function img = load_tracking(s,ind,type,thepath)
% LOAD_TRACKING Loads images for the tracking assignment.
%
% img = load_tracking(s,i)
% img = load_tracking(s,i,type)
% img = load_tracking(s,i,type,thepath)
%
% s: Which of the three sequences to load.
% i: Vector of image(s) indeces, possible ranges are:
%    s=1: 1-79
%    s=2: 1-98
%    s=3: 1-69
% type: 'p' to load an image containing foreground probabilities
%       (default)
%       'i' to load the original (RGB) image 
% thepath: The path to the *.png-files.  By default, we assume
%          that the files are located in the same directory as when
%          installed on the EFD computer network, i.e.
%          ./data/protected/streets/torna
%          If the fucntion can't find an image file it will revert to the
%          default path for a second attempt before failing.
%
% img: The image(s), converted to double and scaled to 0-1.
%      Multiple images are return as a 4-D array for RGB images and 3-D
%      array for foreground probabilities.
%
% Short description of the sequences.
% Sequence 1: Car entering at the lower right corner and making a
%             left turn.
% Sequence 2: Several cars travelling right to left with numerous
%             overlapping foreground estimates. Especially the green
%             car making a left turn would be interesting to
%             track.
% Sequence 3: Two cars entering from the right, one heading across
%             the intersection and one making a right turn. 

% Copyright (c) 2005 Johan Lindström
% $Id: load_tracking.m 3006 2006-10-05 14:25:17Z johanl $

if(nargin<3), type = []; end;
if(nargin<4), thepath = []; end;

if(isempty(type))
  type = 'p';
end;

if (isempty(thepath))
  [p,n,e]=fileparts(which(mfilename));
  thepath = fullfile(p,filesep,'data',filesep,'protected',filesep, ...
		     'streets',filesep,'torna',filesep); 
end

if(type=='p')
  filename = sprintf('Pbg%d_',s);
  img = zeros(251,251,length(ind));
elseif(type=='i')
  filename = sprintf('img%d_',s);
  img = zeros(251,251,3,length(ind));
else
  error( sprintf('Unknown type: %s',type) );
end;

for k=1:length(ind)
  try
    tmp = imread(sprintf('%s%s%03d.png',thepath,filename,ind(k)));
  catch %failed, try default path.
    [p,n,e]=fileparts(which(mfilename));
    thepath = fullfile(p,filesep,'data',filesep,'protected',filesep, ...
		       'streets',filesep,'torna',filesep); 
    tmp = imread(sprintf('%s%s%03d.png',thepath,filename,ind(k)));
  end;
  if(type=='p')
    img(:,:,k) = double(tmp)/255;
  else
    img(:,:,:,k) = double(tmp)/255;
  end;
end;
