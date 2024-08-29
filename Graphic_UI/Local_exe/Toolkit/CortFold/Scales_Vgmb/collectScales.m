function SubjectDataTable=collectScales(targetpath,fpath,scales,lrstr)

%this function basically assembles an overview table for the subject across
%all scales. It also adds a scale 0 for the original mesh (where available)

%INPUTS:

%-targetpath is a string of the path where the outputs from estimateScale 
% were saved. Further outputs will also be stored here.

%-fpath is the path of the Freesurfer surf folder of the subject.

%-scales is a vector of all the scales you want to assemble into one table

%-lrstr is a string of 'l' or 'r' to indicate the hemisphere to run





outputname=[targetpath '/AllScales_hemi=' lrstr '.mat'];
xlsname=[targetpath '/AllScales_hemi=' lrstr '.xls'];


%% make table
z=zeros(size(scales));
At=z;Ae=z;CH=z;GMV=z;nTriAt=z;WMarea=z;At_woCC=z;CH_woCC=z;
for k=1:length(scales)
    scale=scales(k)
    voxname=[targetpath '/VoxelisationVol_scale=' num2str(scale) '_hemi=' lrstr '.mat'];
    
    if exist([voxname])==2 
        load(voxname)
        %could read out more of the output, but this is enough for now
        nTriAt(k)=output.mesh_equinTri_pialFull;
        At(k)=output.mesh_area_pialFull;
        CH(k)=output.mesh_area_CHFull;
        GMV(k)=output.volume_GM;
        WMarea(k)=output.mesh_area_whiteFull;
        if isfield(output,'mesh_area_pial')
        At_woCC(k)=output.mesh_area_pial;
        CH_woCC(k)=output.mesh_area_CH;
        end

    else
        nTriAt(k)=NaN;
        At(k)=NaN;
        CH(k)=NaN;
        GMV(k)=NaN;
        WMarea(k)=NaN;
        warning('Voxelisation file does not exist! Filled with NaNs')
    end
    

end

if isfield(output,'mesh_area_pial')
    SubjectDataTable=table(scales',GMV',nTriAt',At',WMarea',CH',Ae',At_woCC',CH_woCC');
else
    SubjectDataTable=table(scales',GMV',nTriAt',At',WMarea',CH',Ae');
end
SubjectDataTable.Properties.VariableNames{1}='Scale';
SubjectDataTable.Properties.VariableNames{2}='GM_Vol';
SubjectDataTable.Properties.VariableNames{3}='n_Tri';
SubjectDataTable.Properties.VariableNames{4}='At';
SubjectDataTable.Properties.VariableNames{5}='WM_area';
SubjectDataTable.Properties.VariableNames{6}='CH';
SubjectDataTable.Properties.VariableNames{7}='Ae';%Ae is all zeros, as we don't calculate it currently. Could do in future. Keeping it for now as a placeholder. At scale=0 we do have a value for Ae if FS folder available.
if isfield(output,'mesh_area_pial')
    SubjectDataTable.Properties.VariableNames{8}='At_woCC';
    SubjectDataTable.Properties.VariableNames{9}='CH_woCC';
end



%if FS folder is provided
if ~isempty(fpath) && exist([fpath lrstr 'h.pial'],'file')==2 && exist([fpath lrstr 'h.white'],'file')==2 && exist([fpath lrstr 'h.pial-outer-smoothed'],'file')==2

    SubjectDataTable{length(scales)+1,1}=0;

    [pialv,pialf]=freesurfer_read_surf([fpath lrstr 'h.pial']);
    [whitev,whitef]=freesurfer_read_surf([fpath lrstr 'h.white']);
    [opialv,opialf]=freesurfer_read_surf([fpath lrstr 'h.pial-outer-smoothed']); 
    volp=calcMeshVol(pialf,pialv);
    volw=calcMeshVol(whitef,whitev);
    SubjectDataTable{length(scales)+1,2}=volp-volw;

    SubjectDataTable{length(scales)+1,3}=Inf;

    At_orig=calcTriangleArea(pialf,pialv);
    WM_orig=calcTriangleArea(whitef,whitev);
    Ae_orig=calcTriangleArea(opialf,opialv);
    K=convhull(pialv);
    CH_orig=calcTriangleArea(K,pialv);
    SubjectDataTable{length(scales)+1,4}=sum(At_orig);
    SubjectDataTable{length(scales)+1,5}=sum(WM_orig);
    SubjectDataTable{length(scales)+1,6}=sum(CH_orig);
    SubjectDataTable{length(scales)+1,7}=sum(Ae_orig);
else
    SubjectDataTable{length(scales)+1,1}=0;
    SubjectDataTable{length(scales)+1,2:7}=NaN;
    if isfield(output,'mesh_area_pial')
        SubjectDataTable{length(scales)+1,2:9}=NaN;
    end

    warning('Freesurfer folder does not exist, filled scale=0 with NaNs')
end

writetable(SubjectDataTable,xlsname);
save(outputname,'SubjectDataTable');

