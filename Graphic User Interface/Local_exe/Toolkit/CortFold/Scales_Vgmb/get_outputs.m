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
    
    %for each voxel centre find the nearest surface labelnumNaNs
    vid=find(M>0);
    [vid_x,vid_y,vid_z]=ind2sub(size(M),vid);
    pvv=[vid_y,vid_x,vid_z]*scale+[xmin ymin zmin]-scale/2;%convert to the same space as pialv, don't ask me why x and y are swapped, the indexing function swaps it.
    pvv_label=matchSurfLabel(cclabel,pialv,pvv);
    rmids=find(pvv_label==0);
    ind = sub2ind(size(M),vid_x(rmids),vid_y(rmids),vid_z(rmids));
    M(ind)=0;

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f_pial,v_pial]=isosurface(in_pial,.5);

output = table();

output.Scale = scale;

output.AvgCortThickness  = nan;
output.PialArea          = nan;
output.SmoothPialArea    = nan;
output.WhiteArea         = nan;
output.GaussianCurvature = nan;
output.ConvexHullArea    = nan;

output.PialFullArea      = nan;
output.WhiteFullArea     = nan;
output.SmoothPialFullArea= nan;
output.ConvexHullFullArea= nan;

output.PialVol           = nan;
output.WhiteVol          = nan;
output.GreymatterVol     = nan;

output.PialFullVol       = nan;
output.WhiteFullVol      = nan;
output.GreymatterFullVol = nan;
output.SmoothPialFullVol = nan;

if ~isempty(v_pial)

    %shift back to original coordinates
    v_pial=v_pial*scale+[xmin ymin zmin]-scale/2;

    if cce==1
        label_mesh_pial=matchSurfLabel(cclabel,pialv,v_pial);
        cca=calcPartArea(f_pial,v_pial,find(label_mesh_pial==0));
        f_pial_only=discardCCsubc(f_pial,label_mesh_pial);
    end

    area=calcTriangleArea(f_pial,v_pial);
    output.PialFullArea=sum(area);
    
    %CH of this iteration
    K=convhull(v_pial);
    area=calcTriangleArea(K,v_pial);
    output.ConvexHullFullArea=sum(area);

    if cce==1 %if cclabel data is available
        output.PialArea=output.PialFullArea-cca;
        hull_curv = ismember(K, find(label_mesh_pial==1));
        K_only = K(min(hull_curv, [], 2) == 1,:);
        output.ConvexHullArea = sum(calcTriangleArea(K_only, v_pial));
    end

    [f_white,v_white]=isosurface(in_white,.01);%as all the voxels are fully inside the WM (unlike the GM surface)

    if ~isempty(f_white)
        v_white=v_white*scale+[xmin ymin zmin]-scale/2;

        area=calcTriangleArea(f_white,v_white);
        output.WhiteFullArea=sum(area);

        if cce==1
            output.WhiteArea=output.WhiteFullArea-cca;
        end
    end    
end
end