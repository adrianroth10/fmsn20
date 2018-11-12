function fms150path
% FMS150PATH Set path to fms150-subdirectories.
%
% Is called automatically by fms150 on the efd network.
% Call manually on other systems.

% $Id: fms150path.m 2989 2006-09-28 07:21:56Z johanl $

[p,n,e]=fileparts(which(mfilename));
addpath( genpath([p]) );
% $$$ addpath([p]);
% $$$ addpath([p,filesep,'data']);
% $$$ addpath([p,filesep,'data',filesep,'markov']);
% $$$ addpath([p,filesep,'data',filesep,'face']);
% $$$ addpath([p,filesep,'data',filesep,'face',filesep,'2004']);
% $$$ addpath([p,filesep,'data',filesep,'puzzle']);
% $$$ addpath([p,filesep,'data',filesep,'face',filesep,'2005']);
% $$$ addpath([p,filesep,'data',filesep,'face',filesep,'bioid']);
% $$$ addpath([p,filesep,'data',filesep,'photo']);
% $$$ addpath([p,filesep,'data',filesep,'streets',filesep,'torna']);
% $$$ addpath([p,filesep,'classification']);
% $$$ addpath([p,filesep,'cov']);
% $$$ addpath([p,filesep,'graphics']);
% $$$ addpath([p,filesep,'markov']);
% $$$ addpath([p,filesep,'misc']);
% $$$ addpath([p,filesep,'shape']);
% $$$ addpath([p,filesep,'warp']);

set(0,'defaultfiguredoublebuffer','on')
