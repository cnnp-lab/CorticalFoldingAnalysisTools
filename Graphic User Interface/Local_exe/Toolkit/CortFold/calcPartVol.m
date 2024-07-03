function TotalVol=calcPartVol(face,vertex,vertexid)

[liaa]=ismember(face,vertexid);
fid=sum(liaa,2);
    
vol3=calcMeshVol(face(fid==3,:),vertex);
vol2=calcMeshVol(face(fid==2,:),vertex);
vol1=calcMeshVol(face(fid==1,:),vertex);
    
TotalVol=sum(vol3)+sum(vol2)*2/3+sum(vol1)*1/3;