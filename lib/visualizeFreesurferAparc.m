function visualizeFreesurferAparc(figNumber,stitle,valuePerROI,refNames,cllim,mycmap_hotcold)

% valuePerROI is a vector of 68 long, left hemisphere first (1-34)
% refNames is a list of names (categorical) of the ROIs, left first (1-34)


load fsavg_surf_label.mat %making this file has gone into ~/Syncthing/FoldingScaling/scripts/make_fsavg_surf_label

%left first!
ROIvalue=valuePerROI(1:34);
ROInames=refNames(1:34);

%match ROI names left
ROIcolorcode=zeros(length(ROInames),1);
for k=1:length(lct.struct_names)
    mid=find(ROInames==lct.struct_names{k});
    ROIcolorcode(mid)=lct.table(k,5);
end

%make colorvector
ROI2llabel_ds=zeros(length(llabel_ds),1);
for k=1:length(ROIcolorcode)
    rc=ROIcolorcode(k);
    ROI2llabel_ds(llabel_ds==rc)=ROIvalue(k);
end


%right second
ROIvalue=valuePerROI(35:end);
ROInames=refNames(35:end);

%match ROI names left
ROIcolorcode=zeros(length(ROInames),1);
for k=1:length(rct.struct_names)
    mid=find(ROInames==rct.struct_names{k});
    ROIcolorcode(mid)=rct.table(k,5);
end

%make colorvector
ROI2rlabel_ds=zeros(length(rlabel_ds),1);
for k=1:length(ROIcolorcode)
    rc=ROIcolorcode(k);
    ROI2rlabel_ds(rlabel_ds==rc)=ROIvalue(k);
end




figure(figNumber)
subplot(2,4,1)
trisurf(lpialf_ds,lpialv_ds(:,1),lpialv_ds(:,2),lpialv_ds(:,3),ROI2llabel_ds,'EdgeAlpha',0.2)
caxis(cllim)
colormap(mycmap_hotcold)
view(-90, 0)
axis tight
box off
axis off
set(gca,'color','none') 
title('Left')

subplot(2,4,5)
trisurf(lpialf_ds,lpialv_ds(:,1),lpialv_ds(:,2),lpialv_ds(:,3),ROI2llabel_ds,'EdgeAlpha',0.2)
caxis(cllim)
colormap(mycmap_hotcold)
view(90, 0)
axis tight
box off
axis off
set(gca,'color','none') 

subplot(2,4,[2 3 6 7])
trisurf(lpialf_ds,lpialv_ds(:,1),lpialv_ds(:,2),lpialv_ds(:,3),ROI2llabel_ds,'EdgeAlpha',0.2)
hold on
trisurf(rpialf_ds,rpialv_ds(:,1),rpialv_ds(:,2),rpialv_ds(:,3),ROI2rlabel_ds,'EdgeAlpha',0.2)
hold off
caxis(cllim)
colormap(mycmap_hotcold)
view(0, 90)
box off
axis off
set(gca,'color','none') 

subplot(2,4,4)
trisurf(rpialf_ds,rpialv_ds(:,1),rpialv_ds(:,2),rpialv_ds(:,3),ROI2rlabel_ds,'EdgeAlpha',0.2)
caxis(cllim)
colormap(mycmap_hotcold)
view(90, 0)
axis tight
box off
axis off
set(gca,'color','none') 
title('Right')

subplot(2,4,8)
trisurf(rpialf_ds,rpialv_ds(:,1),rpialv_ds(:,2),rpialv_ds(:,3),ROI2rlabel_ds,'EdgeAlpha',0.2)
caxis(cllim)
colormap(mycmap_hotcold)
view(-90, 0)
axis tight
box off
axis off

suptitle(stitle)
set(gca,'color','none') 
