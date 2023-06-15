function equi_n=recalc_nTri_uniformTriangles(scale,area)

tris=sqrt(2*(scale/2)^2)*scale/2;
ua=unique(area);
n=zeros(length(ua),1);
for k=1:length(ua)
    fid=find(area==ua(k));
    n(k)=numel(fid)*ua(k)/tris;
    
end

equi_n=sum(n);