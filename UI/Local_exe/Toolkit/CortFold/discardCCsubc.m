function f_pial_only=discardCCsubc(f_pial,label_mesh_pial)

keepid=find(label_mesh_pial==1);
[liaa]=ismember(f_pial,keepid);
f_pial_only=f_pial;
f_pial_only(sum(liaa,2)<3,:)=[];