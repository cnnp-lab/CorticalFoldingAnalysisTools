function MSfb=makeFacebased(MSvb,faces)
    MSfb=zeros(size(faces,1),1);
    for k=1:size(faces,1)
        threeValues=MSvb(faces(k,:));
        if isempty(find(threeValues==0,1))
            MSfb(k)=mean(MSvb(faces(k,:)));
        else
            MSfb(k)=0;
        end
    end
end