function output=get_outputs(in_pial, in_white, xmin, ymin, zmin, scale, pialv,cclabel)
M=(in_pial-in_white);
M(M<0)=0;

if isempty(cclabel)
    cce=0;
else
    cce=1;
end

%remove any voxels in M on the cc (if info is available)
%to do, M is currently not saved, should really save it!
if cce==1
disp('removing voxels around corpus collosum')
tic

    %for each voxel centre find the nearest surface labelnumNaNs
    vid=find(M>0);
    [vid_x,vid_y,vid_z]=ind2sub(size(M),vid);
    pvv=[vid_y,vid_x,vid_z]*scale+[xmin ymin zmin]-scale/2;%convert to the same space as pialv, don't ask me why x and y are swapped, the indexing function swaps it.
    pvv_label=matchSurfLabel(cclabel,pialv,pvv);    
    rmids=find(pvv_label==0);
    ind = sub2ind(size(M),vid_x(rmids),vid_y(rmids),vid_z(rmids));
    M(ind)=0;
    toc
end
%         %test out code above works:
%         figure()
%         trisurf(whitef,whitev(:,1),whitev(:,2),whitev(:,3),'FaceColor',[0.8 0.8 0.8])
%         hold on
%         %trisurf(f_pial,v_pial(:,1),v_pial(:,2),v_pial(:,3))
%         pvv=[vid_y,vid_x,vid_z]*scale+[xmin ymin zmin]-scale/2;%convert to the same space as pialv
%         scatter3(pvv(:,1),pvv(:,2),pvv(:,3),30,pvv_label,'filled')
%         hold off
%         caxis([0 1])
%         

output.volume_GM=sum(sum(sum(M)))*scale^3;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('get At surface area and cc area')
tic
[f_pial,v_pial]=isosurface(in_pial,.5);

if ~isempty(v_pial)

%shift back to original coordinates
v_pial=v_pial*scale+[xmin ymin zmin]-scale/2;

if cce==1
    label_mesh_pial=matchSurfLabel(cclabel,pialv,v_pial);
    cca=calcPartArea(f_pial,v_pial,find(label_mesh_pial==0));
    f_pial_only=discardCCsubc(f_pial,label_mesh_pial);
end

area=calcTriangleArea(f_pial,v_pial);
output.mesh_area_pialFull=sum(area);
output.mesh_nTri_pialFull=size(f_pial,1);
output.mesh_avgSizeTri_pialFull=mean(area);
%output.mesh_equinTri_pialFull=recalc_nTri_uniformTriangles(scale,area);
output.mesh_equinTri_pialFull=NaN;%not used
toc


disp('get Ae or CH surface area')
tic
%CH of this iteration
K=convhull(v_pial);
area=calcTriangleArea(K,v_pial);
output.mesh_area_CHFull=sum(area);
toc

disp('match cc to CH surface area')
tic
if cce==1 %if cclabel data is available
    output.mesh_area_pial=output.mesh_area_pialFull-cca;
    hull_curv = ismember(K, find(label_mesh_pial==1));
    K_only = K(min(hull_curv, [], 2) == 1,:);
    output.mesh_area_CH = sum(calcTriangleArea(K_only, v_pial));

    area=calcTriangleArea(f_pial_only,v_pial);
    output.mesh_avgSizeTri_pial=mean(area);
    %output.mesh_equinTri_pial=recalc_nTri_uniformTriangles(scale,area);
    output.mesh_nTri_pial=size(f_pial_only,1);
else
    output.mesh_area_pial=NaN;
    output.mesh_area_CH=NaN;

    output.mesh_avgSizeTri_pial=NaN;
    %output.mesh_equinTri_pial=NaN;
    output.mesh_nTri_pial=NaN;
end
output.mesh_equinTri_pial=NaN;%not used

toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Get WM surface area')
tic
[f_white,v_white]=isosurface(in_white,.01);%as all the voxels are fully inside the WM (unlike the GM surface)
%[v_white,f_white]=binsurface(in_white,3);%iso2mesh equivalent of the above
if ~isempty(f_white)
    v_white=v_white*scale+[xmin ymin zmin]-scale/2;

    area=calcTriangleArea(f_white,v_white);
    output.mesh_area_whiteFull=sum(area);

    output.mesh_nTri_whiteFull=size(f_white,1);
    output.mesh_avgSizeTri_whiteFull=mean(area);
    output.mesh_equinTri_whiteFull=recalc_nTri_uniformTriangles(scale,area);
    
    if cce==1
        output.mesh_area_white=output.mesh_area_whiteFull-cca;
    else
        output.mesh_area_white=NaN;
    end
    
else
    output.mesh_area_whiteFull=NaN;

    output.mesh_nTri_whiteFull=NaN;
    output.mesh_avgSizeTri_whiteFull=NaN;
    output.mesh_equinTri_whiteFull=NaN;
end
toc

else
    output.mesh_area_pialFull=NaN;
    output.mesh_nTri_pialFull=NaN;
    output.mesh_avgSizeTri_pialFull=NaN;
    output.mesh_equinTri_pialFull=NaN;%not used
    output.mesh_area_pial=NaN;
    output.mesh_area_CH=NaN;
    output.mesh_area_CHFull=NaN;
    output.mesh_avgSizeTri_pial=NaN;
    %output.mesh_equinTri_pial=NaN;
    output.mesh_nTri_pial=NaN;
    output.mesh_equinTri_pial=NaN;%not used
    output.mesh_area_white=NaN;
    output.mesh_area_whiteFull=NaN;
    output.mesh_nTri_whiteFull=NaN;
    output.mesh_avgSizeTri_whiteFull=NaN;
    output.mesh_equinTri_whiteFull=NaN;
end
end