clear all
clc

% Apply all the scales and iterations
targetScale=[0.225:0.025:.375];%[0.225:0.025:0.3000]
nIter=[5 5 4 4 4 4 4];%[ ];

lrstr='lr';

cce = 0;

hcpsubjs='0010001';
root = '/Users/guillermo/Sync/Data/04_CorticalFolding_Interface/gmb_Trial/FreeSurfer/';

addpath(genpath('/Applications/freesurfer/7.4.1'));

subjpath=fullfile(root, hcpsubjs);
scaled_Outputs = fastEstimateScale(targetScale,nIter,subjpath,lrstr,cce);
