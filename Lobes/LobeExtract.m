function LobeExtract(dataset,rootfolder,subjrootfolder,sl,POSname,targetpath)

lrstr='lr';

if strcmp(dataset,'ABIDE2')
    subfstr='/session_1/FS_OUTPUT/';%just for ABIDE2
else
    
    subfstr='';
end

addpath([rootfolder 'lib/'])%the lib files are attached
addpath([rootfolder 'lib/FSmatlab/'])
addpath([rootfolder 'lib/iso2mesh/'])

load([rootfolder 'data_collated/' dataset '.mat'])
load(['LUT_lobes.mat'])

%Note: labels 0 1 2 3 4 5 are CC F P T O Insula

%%
eval(['subjid=' dataset '.SubjID;']);


z=zeros(6,2);
AvgThickness=z;
TotalArea=z;
SmoothArea=z;
K_ds=z;
SmoothArea_ds=z;
K_CHwoB=z;
SmoothArea_CHwoB=z;

for lr=1:2

    ID=subjid{sl};
    ID=finFSdsubjID(ID,subjrootfolder);

    disp([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) 'h.pial-outer-smoothed'])
    if exist([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) 'h.pial-outer-smoothed'])==2
        [Thickness, ~] = read_curv([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) 'h.thickness']);
        [pialv,pialf]=freesurfer_read_surf([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) 'h.pial']);
        [opialv,opialf]=freesurfer_read_surf([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) POSname]);
        [whitev,whitef]=freesurfer_read_surf([subjrootfolder '/' ID subfstr '/surf/' lrstr(lr) 'h.white']);
        [~,label,~]=read_annotation([subjrootfolder '/' ID subfstr '/label/' lrstr(lr) 'h.aparc.annot']);
        
        %% match labels to smooth surface

        tic
        label_smooth=matchSurfLabel(label,pialv,opialv);
        toc
        
        %% relabel for lobes
        tic
        nlabel_total_lobe=label;
        nlabel_smooth_lobe=label_smooth;
        for k=0:5
            idtoget=LUT_lobes(LUT_lobes(:,2)==k,1);

            for l=1:length(idtoget)
                nlabel_total_lobe(nlabel_total_lobe==idtoget(l))=k;
                nlabel_smooth_lobe(nlabel_smooth_lobe==idtoget(l))=k;
            end

        end
        toc
        
        %% extract measures for lobes
        
        tic
        [AvgThicknessLR,TotalAreaLR,SmoothAreaLR,~]=getBasicMeasuresForParts(nlabel_total_lobe,nlabel_smooth_lobe,pialf,pialv,opialf,opialv,Thickness,whitef,whitev);
        toc
        
        %% calculate K from downsampled opial 
        tic
        %downsample Ae
        chr=0.2;
        [opialv_ds,opialf_ds]=meshresample(opialv,opialf,chr);
        
        %reassign labels
        nlabel_smooth_lobe_downsampled=matchSurfLabel(nlabel_smooth_lobe,opialv,opialv_ds);
        
        %calculate area
        [dsAreaCapi,~]=calcPartAreai(opialf_ds,opialv_ds,1:size(opialv_ds,1));
        
        %get curvature
        Gnode=getGaussianCurvPart(opialf_ds,opialv_ds,1:length(opialv_ds));

        KLR_frm_ds_opial=zeros(6,1);
        Area_frm_ds_opial=zeros(6,1);
        for l=0:5
            KLR_frm_ds_opial(l+1)=sum(Gnode(nlabel_smooth_lobe_downsampled==l));
            Area_frm_ds_opial(l+1)=sum(dsAreaCapi(nlabel_smooth_lobe_downsampled==l));
        end

        toc
        
        %% use CH of opial_ds segments & ignore boundary
        
        tic
        KLR_frm_CHof_dsopial=zeros(6,1);
        Area_frm_CHof_dsopial=zeros(6,1);
        for l=0:5
            ovids=find(nlabel_smooth_lobe_downsampled==l);
            Lobepoints=opialv_ds(ovids,:);
            if length(Lobepoints)>4
            CHSf_ds = convhull(Lobepoints);
            Gru=getGaussianCurvPart(CHSf_ds,Lobepoints,1:length(ovids));

            isBoundary=zeros(length(ovids),1);
            for kl=1:length(ovids)
                %where does idFL pop up in opialf?
                fid=opialf_ds(:,1)==ovids(kl) | opialf_ds(:,2)==ovids(kl) | opialf_ds(:,3)==ovids(kl);

                neighbourIDs=unique(opialf_ds(fid,:));
                ncolors=numel(unique(nlabel_smooth_lobe_downsampled(neighbourIDs)));
                if ncolors>1
                    isBoundary(kl)=1;
                end

            end
            
            GrCHS=sum(Gru(isBoundary==0));%read out curv
            Bids=find(isBoundary==1);
            [TotalAreaCapi,~]=calcPartAreai(CHSf_ds,Lobepoints,1:size(Lobepoints,1));
            Btri=ismember(CHSf_ds(:,1),Bids) & ismember(CHSf_ds(:,2),Bids) & ismember(CHSf_ds(:,3),Bids);
            CHSArea=sum(TotalAreaCapi(~Btri));%read out surface again
            
            KLR_frm_CHof_dsopial(l+1)=GrCHS;
            Area_frm_CHof_dsopial(l+1)=CHSArea;
            end
        end
        toc
        
        %% write out
        AvgThickness(:,lr)=AvgThicknessLR;
        TotalArea(:,lr)=TotalAreaLR;
        SmoothArea(:,lr)=SmoothAreaLR;
        K_ds(:,lr)=KLR_frm_ds_opial;
        SmoothArea_ds(:,lr)=Area_frm_ds_opial;
        K_CHwoB(:,lr)=KLR_frm_CHof_dsopial;
        SmoothArea_CHwoB(:,lr)=Area_frm_CHof_dsopial;
        
    end
    
end


save([targetpath dataset '_LobesExtract' num2str(sl) '.mat'],'AvgThickness','TotalArea','SmoothArea','K_ds','K_CHwoB','SmoothArea_ds','SmoothArea_CHwoB','sl')

