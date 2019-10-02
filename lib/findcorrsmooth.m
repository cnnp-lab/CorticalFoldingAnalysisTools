function idssmooth=findcorrsmooth(faceids,smoothfaces,smoothvertices,faces,vertices)
smoothT=triangulation(double(smoothfaces),double(smoothvertices));
normalT=triangulation(double(faces),double(vertices));

smoothIC = incenter(smoothT);
normalIC = incenter(normalT);

allmind=zeros(length(faceids),3);
parfor k=1:length(faceids)
    d=sqrt((smoothIC(:,1)-normalIC(faceids(k),1)).^2 + (smoothIC(:,2)-normalIC(faceids(k),2)).^2 + (smoothIC(:,3)-normalIC(faceids(k),3)).^2);
    [v,sid]=sort(d);
    if v(1)<5;
        allmind(k,:)=sid(1:3);
    end

end

idssmooth=unique(allmind);
idssmooth(idssmooth==0)=[];
