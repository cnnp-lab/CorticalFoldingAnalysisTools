function fid=vid2fid(vids,faces)

c1=ismember(faces(:,1),vids);
c2=ismember(faces(:,2),vids);
c3=ismember(faces(:,3),vids);

fid=find(c1 & c2 & c3);

end