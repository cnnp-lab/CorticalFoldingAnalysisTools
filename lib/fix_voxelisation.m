function in_pial_fixed=fix_voxelisation(in_pial)

%% kernels for detecting rods
intkernels=[];
intkernels(:,:,1)=[-1 -1 -1;
                  -1 -1 -1;
                  -1 -1 -1];
intkernels(:,:,2)=[1 1 1;
                  1 1 1;
                  1 1 1];

              
in_pial_cvds=convn(in_pial,intkernels,'same');

%% FIX bottom
[fidx,fidy]=find(in_pial(:,:,2)==1 & in_pial(:,:,3)==1); %potential start of rod

in_pial_fixed=in_pial;

for k=1:length(fidx)
%     figure(1)
%     plot(squeeze(in_pial(fidx(k),fidy(k),:)),'--')
%     hold on
%     plot(squeeze(in_pial_cvd(fidx(k),fidy(k),:)))
%     plot(squeeze(in_pial_cvds(fidx(k),fidy(k),:)))
    
    
    stop_rod_id=min(find(squeeze(in_pial_cvds(fidx(k),fidy(k),:))<0));
    
    in_pial_fixed(fidx(k),fidy(k),1:stop_rod_id-1)=0;
    
%     plot(squeeze(in_pial_fixed(fidx(k),fidy(k),:)),'.')
%     hold off
%     pause()
end


%% FIX top
[fidx,fidy]=find(in_pial(:,:,end-2)==1 & in_pial(:,:,end-1)==1); %potential start of rod



for k=1:length(fidx)
%     figure(1)
%     plot(squeeze(in_pial(fidx(k),fidy(k),:)),'--')
%     hold on
%     plot(squeeze(in_pial_cvd(fidx(k),fidy(k),:)))
%     plot(squeeze(in_pial_cvds(fidx(k),fidy(k),:)))
    
    
    stop_rod_id=find(squeeze(in_pial_cvds(fidx(k),fidy(k),:))>0);
    stop_rod_id(stop_rod_id==size(in_pial,3))=[];
    stop_rod_id(stop_rod_id==size(in_pial,3)-1)=[];
    stop_rod_id=max(stop_rod_id);
    
    in_pial_fixed(fidx(k),fidy(k),stop_rod_id:end)=0;
    
%     plot(squeeze(in_pial_fixed(fidx(k),fidy(k),:)),'.')
%     hold off
%     pause()
end

