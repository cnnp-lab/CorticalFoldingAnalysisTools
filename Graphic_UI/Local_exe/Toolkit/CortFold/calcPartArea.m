function TotalArea=calcPartArea(face,vertex,vertexid)

[liaa]=ismember(face,vertexid);
fid=sum(liaa,2);
    
area3=calcTriangleArea(face(fid==3,:),vertex);
area2=calcTriangleArea(face(fid==2,:),vertex);
area1=calcTriangleArea(face(fid==1,:),vertex);
    
TotalArea=sum(area3)+sum(area2)*2/3+sum(area1)*1/3;