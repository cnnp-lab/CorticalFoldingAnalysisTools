function vol=calcMeshVol(f,v)
    
vol=0;


    for k=1:length(f)
        fk=f(k,:);
        v1=v(fk(1),:);
        v2=v(fk(2),:);
        v3=v(fk(3),:);
        vk=1/6*dot(cross(v1,v2),v3);
        vol=vol+vk;
    end

end