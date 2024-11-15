function [tbl,corrupt] = fastEstimateScale(targetScale,nIter,path,hemis_md,cce)
% To use with a target scale and multiples of the target scale!
%suggest using:
% targetScales=[0.225:0.025:.375];
% nIters=[5 5 4 4 4 4 4];
%
%INPUTS:
%-targetScale is a scalar number indicating the spatial scale in mm to downsample
%to.
%-nIter is the number of iterations to scale up by,
% startScale = targetScale*2^nIter
%
%-fpath is the path of the Freesurfer surf folder of the subject.
%alternatively the folder that contains a ?h.mat file, only used for NHP.
%
%-targetpath is a string of the path where the outputs should be saved.
%Outputs are the downsampled stats and also the 3D image.
%
%-lrstr is a string of 'l' or 'r' to indicate the hemisphere to run
%
%-overwriteold is 1 or 0 to indicate if an older run should be overwritten
%in targetpath
%
%-cce can be set to zero by default. It determines if we remove the corpus
%collosum from volume and surface calculations. If cce=1, CC is removed
%from all calculations, but this slows down the algorithm immensely (1h to 10h).

switch(hemis_md)
    case "left"
        hemis = 'l';
        leftright = "left";
    case "right"
        hemis = 'r';
        leftright = "right";
    otherwise % avg, sum, both
        hemis = 'lr';
        leftright = ["left" "right"];
end

fpath = fullfile(path,'/surf/');

%% INITIAL SCALE
% load surfaces & make checks
tbl = table();

corrupt = zeros(1,length(hemis));

for i=1:length(hemis)
    lrstr = hemis(i);
    if exist([fpath lrstr 'h.mat'],'file')==2
        load([fpath lrstr 'h.mat'])
    elseif exist([fpath lrstr 'h.pial'],'file')==2 && exist([fpath lrstr 'h.white'],'file')==2
        [pialv,pialf]=freesurfer_read_surf([fpath lrstr 'h.pial']);
        [whitev,whitef]=freesurfer_read_surf([fpath lrstr 'h.white']);
    end


    if exist([fpath lrstr 'h.thickness'],'file')==2 && cce==1
        [Thickness, ~] = read_curv([fpath lrstr 'h.thickness']);
        cclabel=single(Thickness>0);
    else
        cclabel=[];
    end

    if exist('pialf','var')~=1 || exist('pialv','var')~=1
        corrupt(i) = 1;
        continue
    end

    if exist('whitef','var')~=1 || exist('whitev','var')~=1
        corrupt(i) = 1;
        continue
    end

    for j=1:length(targetScale)
        % setup initial grid & check if in WM/GM
        startScale = targetScale(j)*2^nIter(j);
        bordersizemin=startScale; %in mm
        bordersize=bordersizemin+3*startScale;%so border is at least 4*scale away from mesh

        xmin=floor(min(pialv(:,1)))-bordersize;
        xmax=ceil(max(pialv(:,1)))+bordersize;
        ymin=floor(min(pialv(:,2)))-bordersize;
        ymax=ceil(max(pialv(:,2)))+bordersize;
        zmin=floor(min(pialv(:,3)))-bordersize;
        zmax=ceil(max(pialv(:,3)))+bordersize;

        %setting up voxelisation
        [xg,yg,zg]=meshgrid(single(xmin:startScale:xmax),single(ymin:startScale:ymax),single(zmin:startScale:zmax));

        % Check all points for starting res WM
        in_whitegrid=inpolyhedron(whitef,whitev,[xg(:),yg(:),zg(:)]);
        in_whitegrid=reshape(in_whitegrid,size(xg));

        % Starting res pial
        in_pialgrid = in_whitegrid;
        check = in_pialgrid==0; % Only check the points which are not in WM
        in_pialgrid(check)=inpolyhedron(pialf,pialv,[xg(check),yg(check),zg(check)]);
        clear check

        % poll initial grid in GM & WM
        %WM:
        poll=getPoll(in_whitegrid);
        in_white=(poll==8);
        clear poll
        in_white=fix_voxelisation(in_white);%fixing "rods" in z dimension as an artifact of bad meshes

        %GM:
        poll=getPoll(in_pialgrid);
        in_pial=poll>=4;%at least 50% of corners need to be in surface to count the voxel
        clear poll
        in_pial=fix_voxelisation(in_pial);%fixing "rods" in z dimension as an artifact of bad meshes

        % save initial scale
        tbl_2=get_outputs(in_pial, in_white, xmin, ymin, zmin, startScale, pialv,cclabel);

        %% iterate over next scales
        for kk = 1:nIter(j)

            scale = startScale/(2^kk);

            in_whitegrid_prev=in_whitegrid;
            in_pialgrid_prev = in_pialgrid;

            % Double the mesh
            [xg2,yg2,zg2]=meshgrid(single(xmin:scale:xmax),single(ymin:scale:ymax),single(zmin:scale:zmax));

            % Truncate so end points are the same as for smaller scale
            yg2 = yg2(1:(end - size(xg2,1) + size(xg,1)*2-1), 1:(end - size(xg2,2) + size(xg,2)*2-1), 1:(end - size(xg2,3) + size(xg,3)*2-1));
            zg2 = zg2(1:(end - size(xg2,1) + size(xg,1)*2-1), 1:(end - size(xg2,2) + size(xg,2)*2-1), 1:(end - size(xg2,3) + size(xg,3)*2-1));
            xg2 = xg2(1:(end - size(xg2,1) + size(xg,1)*2-1), 1:(end - size(xg2,2) + size(xg,2)*2-1), 1:(end - size(xg2,3) + size(xg,3)*2-1));

            %wm grid---------------------------------------
            in_whitegrid = false(size(xg2));
            in_whitegrid(1:2:end, 1:2:end, 1:2:end) = in_whitegrid_prev;

            % erode x 3
            se=strel('sphere',3);
            dilateI = imdilate(in_whitegrid,se);

            % backfill core WM
            se=strel('sphere',2);
            dilateIc = imdilate(in_whitegrid,se);
            se=strel('sphere',3);%setting this to two might be enough and quicker!
            erodeI=imerode(dilateIc,se);

            % only recalc eroded voxels again
            check=dilateI+in_whitegrid+erodeI==1;
            in_whitegrid(check)=inpolyhedron(whitef,whitev,[xg2(check),yg2(check),zg2(check)]);
            in_whitegrid(erodeI == 1) = 1;
            clear dilateI erodeI dilateIc check

            %poll
            poll=getPoll(in_whitegrid);
            in_white=poll==8;
            clear poll

            %fixing "rods" in z dimension as an artifact of bad meshes
            in_white=fix_voxelisation(in_white);

            %gm grid--------------------------------------
            in_pialgrid = false(size(xg2));
            in_pialgrid(1:2:end, 1:2:end, 1:2:end) = in_pialgrid_prev;
            in_pialgrid(in_whitegrid==1)=1;

            %erode, take out wm
            se=strel('sphere',3);
            dilateI = imdilate(in_pialgrid,se);


            % only recalc eroded voxels again
            check=dilateI+in_pialgrid==1;
            in_pialgrid(check)=inpolyhedron(pialf,pialv,[xg2(check),yg2(check),zg2(check)]);
            clear dilateI
            clear check

            %poll
            poll=getPoll(in_pialgrid);
            in_pial=poll>=4;
            clear poll

            %fixing "rods" in z dimension as an artifact of bad meshes
            in_pial=fix_voxelisation(in_pial);

            % save initial scale
            % NOT store scale 1 to avoid issuess
            output=get_outputs(in_pial, in_white, xmin, ymin, zmin, scale, pialv, cclabel);
            tbl_2 = [tbl_2;output];

            %reset loop
            xg = xg2;

        end

        tbl_2.Hemisphere = repmat(leftright(i),size(tbl_2,1),1);
        tbl = [tbl;tbl_2];
        
    end
end

tbl.Region = repmat("hemisphere",size(tbl,1),1);

tbl = movevars(tbl,"Region",'Before',1);
tbl = movevars(tbl,"Hemisphere",'Before',1);




