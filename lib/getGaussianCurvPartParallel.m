function G=getGaussianCurvPartParallel(pialf,pialv,ids)

G=zeros(length(ids),1);
parfor l=1:length(ids)
    %find list of faces, loop through them
    [fids,fidp]=find(pialf==ids(l));
    if ~isempty(fids)
        alpha=zeros(length(fids),1);
        for f=1:length(fids)
            %find angle
            vids=pialf(fids(f),:);
            mid=mod([fidp(f) fidp(f)+1 fidp(f)+2],3);%rotate the index such that the vertex in question is first
            mid(mid==0)=3;
            vids=vids(mid);

            p1=pialv(vids(1),:);
            p2=pialv(vids(2),:);
            p3=pialv(vids(3),:);

            v2=p1-p2;
            v3=p1-p3;
            alpha(f)=acos(dot(v2,v3)/(norm(v2)*norm(v3)));

        end
        G(l)=2*pi-sum(alpha);
    end
end