function [TotalAreai,fid]=calcPartAreai(face,vertex,vertexid)

[liaa]=ismember(face,vertexid);
fid=sum(liaa,2);


TotalAreai=zeros(length(fid),1);
TotalAreai(fid==3)=calcTriangleArea(face(fid==3,:),vertex);
TotalAreai(fid==2)=2/3*calcTriangleArea(face(fid==2,:),vertex);
TotalAreai(fid==1)=1/3*calcTriangleArea(face(fid==1,:),vertex);
   