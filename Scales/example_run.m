clear all
clc

targetScales=[0.225:0.025:.375];
nIters=[5 5 4 4 4 4 4];


lrstr='lr';

hcpsubjs={'148840'};

%% run across scales
% 

for hn=1:length(hcpsubjs)
    subjpath=['/data1/yujiang/Github/Folding_scales_dev/data/subjects/HCP/' hcpsubjs{hn} '/'];
    targetpath=['/data1/yujiang/Github/Folding_scales_dev/data/subjects/HCP/' hcpsubjs{hn} '/fastESout/'];   
    for lri=1:2
        lr=lrstr(lri);
        for s=1:length(targetScales)
            targetScale=targetScales(s);
            nIter=nIters(s);

            fastEstimateScale(targetScale,nIter,subjpath,targetpath,lr,1,0)
        end
    end
end

%% collect results across scales

scales=[];
for k=1:length(targetScales)
    scales=[scales targetScales(k)*2.^[0:nIters(k)]];
end

scales=sort(scales);

for hn=1:length(hcpsubjs)
    subjpath=['/data1/yujiang/Github/Folding_scales_dev/data/subjects/HCP/' hcpsubjs{hn} '/'];
    targetpath=['/data1/yujiang/Github/Folding_scales_dev/data/subjects/HCP/' hcpsubjs{hn} '/fastESout/'];   
    for lri=1:2
        lr=lrstr(lri);
        SubjectDataTable=collectScales(targetpath,subjpath,scales,lr);
    end
end
