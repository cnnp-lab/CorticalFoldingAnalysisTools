function estimateScale(scale,fpath,lrstr,targetpath,overwriteold,njitter)

%INPUTS: 
%-scale is a scalar number indicating the spatial scale in mm to downsample
%to. ("Downsampling" describes the procedure of coarse-graining the cortex
%with voxels of a particular size here.)

%-fpath is the path of the Freesurfer surf folder of the subject.
%alternatively the folder that contains a ?h.mat file, only used for NHP.

%-lrstr is a string of 'l' or 'r' to indicate the hemisphere to run

%-targetpath is a string of the path where the outputs should be saved.
%Outputs are the downsampled stats and also the 3D image.

%-overwriteold is 1 or 0 to indicate if an older run should be overwritten
%in targetpath

%-njitter is mostly set to 1, which just one iteration of the downsampling
%procedure. If it is set to another number n>1 then the algorith will run n
%iterations, but each time jittering the 3D grid a little. Mainly useful to
%debug/assess robustness.


%% this function is step 1. 
% First generates the voxelisation (for the ribbon) & gets its volume.
% Saves voxelisation & volume info.
% Then measures the At (GM surface) and CH area & numberTriangles -> save

outputname=[targetpath '/VoxelisationVol_scale=' num2str(scale) '_hemi=' lrstr '.mat'];

if exist([outputname],'file')~=2 || overwriteold==1

    
%% load surfaces & make checks

if exist([fpath lrstr 'h.mat'],'file')==2
    load([fpath lrstr 'h.mat'])
elseif exist([fpath lrstr 'h.pial'],'file')==2 && exist([fpath lrstr 'h.white'],'file')==2
    [pialv,pialf]=freesurfer_read_surf([fpath lrstr 'h.pial']);
    [whitev,whitef]=freesurfer_read_surf([fpath lrstr 'h.white']);
end

if exist('pialf','var')~=1 || exist('pialv','var')~=1
    error('No pial surface could be loaded')
end

if exist('whitef','var')~=1 || exist('whitev','var')~=1
    error('No white surface could be loaded')
end


if exist([fpath lrstr 'h.thickness'],'file')==2
    
    [Thickness, ~] = read_curv([fpath lrstr 'h.thickness']);
    cclabel=single(Thickness>0);
    
elseif exist([fpath lrstr 'h.thickness.mat'],'file')==2
    load([fpath lrstr 'h.thickness.mat'])
    cclabel=single(Thickness>0);
end


cce=0;
if exist('cclabel','var') == 1
    cce=1;
else
    warning('No cclabel data loaded')
end


%% start jitter iterations
mkdir(targetpath)
output=struct();
    
    for j=1:njitter
            
        %if we are jittering the grid
        if njitter==1    
            joffset=0;
        else
            joffset=rand(1)*scale;
            disp(['Jitter iteration ' num2str(j)])
        end
        bordersizemin=scale+joffset; %in mm
        bordersize=bordersizemin+3*scale;%so border is at least 4*scale away from mesh


        %% scan scales

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('setup grid system, fill surfaces, get volume')
        tic


        xmin=floor(min(pialv(:,1)))-bordersize;
        xmax=ceil(max(pialv(:,1)))+bordersize;
        ymin=floor(min(pialv(:,2)))-bordersize;
        ymax=ceil(max(pialv(:,2)))+bordersize;
        zmin=floor(min(pialv(:,3)))-bordersize;
        zmax=ceil(max(pialv(:,3)))+bordersize;

        
        
        %setting up voxelisation
        [xg,yg,zg]=meshgrid(single(xmin:scale:xmax),single(ymin:scale:ymax),single(zmin:scale:zmax));
        in_pial=inpolyhedron(pialf,pialv,[xg(:),yg(:),zg(:)]);
        in_pial=reshape(in_pial,size(xg));
        

        
        %transform such that we get have the grid give us 8 corners of a
        %voxel

        poll=zeros(size(in_pial)-1,'uint8');
        poll=uint8(in_pial(1:end-1,1:end-1,1:end-1))...
            +uint8(in_pial(2:end,1:end-1,1:end-1))...
            +uint8(in_pial(2:end,1:end-1,2:end))...
            +uint8(in_pial(1:end-1,2:end,1:end-1))...
            +uint8(in_pial(1:end-1,2:end,2:end))...
            +uint8(in_pial(1:end-1,1:end-1,2:end))...
            +uint8(in_pial(2:end,2:end,1:end-1))...
            +uint8(in_pial(2:end,2:end,2:end));
        in_pial=poll>=4;%at least 50% of corners need to be in surface to count the voxel
        
        %fixing "rods" in z dimension as an artifact of bad meshes
        in_pial=fix_voxelisation(in_pial);
        
        %TO DO: Estimate AE in a different way? - Leave for future
        
        

        in_white=inpolyhedron(whitef,whitev,[xg(:),yg(:),zg(:)]);
        in_white=reshape(in_white,size(xg));
        poll=zeros(size(in_pial)-1,'uint8');
        poll=uint8(in_white(1:end-1,1:end-1,1:end-1))...
            +uint8(in_white(2:end,1:end-1,1:end-1))...
            +uint8(in_white(2:end,1:end-1,2:end))...
            +uint8(in_white(1:end-1,2:end,1:end-1))...
            +uint8(in_white(1:end-1,2:end,2:end))...
            +uint8(in_white(1:end-1,1:end-1,2:end))...
            +uint8(in_white(2:end,2:end,1:end-1))...
            +uint8(in_white(2:end,2:end,2:end));
        in_white=poll==8;

        %fixing "rods" in z dimension as an artifact of bad meshes
        in_white=fix_voxelisation(in_white);
        
        
        M=(in_pial-in_white);
        M(M<0)=0;
        toc
        %%
        disp('removing voxels around corpus collosum')
        tic
        %remove any voxels in M on the cc (if info is available)
        %to do, M is currently not saved, should really save it!
        if cce==1
            %for each voxel centre find the nearest surface label
            vid=find(M>0);
            [vid_x,vid_y,vid_z]=ind2sub(size(M),vid);
            pvv=[vid_y,vid_x,vid_z]*scale+[xmin ymin zmin]-scale/2;%convert to the same space as pialv, don't ask me why x and y are swapped, the indexing function swaps it.
            pvv_label=matchSurfLabel(cclabel,pialv,pvv);    
            rmids=find(pvv_label==0);
            ind = sub2ind(size(M),vid_x(rmids),vid_y(rmids),vid_z(rmids));
            M(ind)=0;
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
        
        output(j).volume_GM=sum(sum(sum(M)))*scale^3;
        toc


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
        output(j).mesh_area_pialFull=sum(area);
        output(j).mesh_nTri_pialFull=size(f_pial,1);
        output(j).mesh_avgSizeTri_pialFull=mean(area);
        %output(j).mesh_equinTri_pialFull=recalc_nTri_uniformTriangles(scale,area);
        output(j).mesh_equinTri_pialFull=NaN;%not used
        toc
        
        
        disp('get Ae or CH surface area')
        tic
        %CH of this iteration
        K=convhull(v_pial);
        area=calcTriangleArea(K,v_pial);
        output(j).mesh_area_CHFull=sum(area);
        toc
        
        disp('match cc to CH surface area')
        tic
        if cce==1 %if cclabel data is available
            output(j).mesh_area_pial=output(j).mesh_area_pialFull-cca;
            hull_curv = ismember(K, find(label_mesh_pial==1));
            K_only = K(min(hull_curv, [], 2) == 1,:);
            output(j).mesh_area_CH = sum(calcTriangleArea(K_only, v_pial));
    
            area=calcTriangleArea(f_pial_only,v_pial);
            output(j).mesh_avgSizeTri_pial=mean(area);
            %output(j).mesh_equinTri_pial=recalc_nTri_uniformTriangles(scale,area);
            output(j).mesh_nTri_pial=size(f_pial_only,1);
        else
            output(j).mesh_area_pial=NaN;
            output(j).mesh_area_CH=NaN;
    
            output(j).mesh_avgSizeTri_pial=NaN;
            %output(j).mesh_equinTri_pial=NaN;
            output(j).mesh_nTri_pial=NaN;
        end
        output(j).mesh_equinTri_pial=NaN;%not used
        
        toc



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Get WM surface area')
        tic
%         [f_white,v_white]=isosurface(in_white,.01);%as all the voxels are fully inside the WM (unlike the GM surface)
        [v_white,f_white]=binsurface(in_white,3);%as all the voxels are fully inside the WM (unlike the GM surface)
        if ~isempty(f_white)
            v_white=v_white*scale+[xmin ymin zmin]-scale/2;

            area=calcTriangleArea(f_white,v_white);
            output(j).mesh_area_whiteFull=sum(area);

            output(j).mesh_nTri_whiteFull=size(f_white,1);
            output(j).mesh_avgSizeTri_whiteFull=mean(area);
            output(j).mesh_equinTri_whiteFull=recalc_nTri_uniformTriangles(scale,area);
            
            if cce==1
                output(j).mesh_area_white=output(j).mesh_area_whiteFull-cca;
            else
                output(j).mesh_area_white=NaN;
            end
            
        else
            output(j).mesh_area_whiteFull=NaN;

            output(j).mesh_nTri_whiteFull=NaN;
            output(j).mesh_avgSizeTri_whiteFull=NaN;
            output(j).mesh_equinTri_whiteFull=NaN;
        end
        toc

        else
            output(j).mesh_area_pialFull=NaN;
            output(j).mesh_nTri_pialFull=NaN;
            output(j).mesh_avgSizeTri_pialFull=NaN;
            output(j).mesh_equinTri_pialFull=NaN;%not used
            output(j).mesh_area_pial=NaN;
            output(j).mesh_area_CH=NaN;
            output(j).mesh_area_CHFull=NaN;
            output(j).mesh_avgSizeTri_pial=NaN;
            %output(j).mesh_equinTri_pial=NaN;
            output(j).mesh_nTri_pial=NaN;
            output(j).mesh_equinTri_pial=NaN;%not used
            output(j).mesh_area_white=NaN;
            output(j).mesh_area_whiteFull=NaN;
            output(j).mesh_nTri_whiteFull=NaN;
            output(j).mesh_avgSizeTri_whiteFull=NaN;
            output(j).mesh_equinTri_whiteFull=NaN;
        end

    end


    save(outputname,'in_pial','in_white','output','scale','njitter','xmin','ymin','zmin');

else
    warning('File already exists and is not overwritten!')
end