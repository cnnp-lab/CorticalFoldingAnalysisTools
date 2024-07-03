function area=calcTriangleArea(face,vertex)

area=zeros(size(face,1),1);
for i=1:size(face,1)
    l=zeros(3,1);
    P1=[vertex(face(i,1),:)];
    P2=[vertex(face(i,2),:)];
    P3=[vertex(face(i,3),:)];
    area(i) = 1/2*norm(cross(P2-P1,P3-P1));

% 
% point1=[vertex(face(i,1),:)];
%     point2=[vertex(face(i,2),:)];
%     point3=[vertex(face(i,3),:)];
%     l(1)=sqrt( (point1(1)-point2(1)).^2 + (point1(2)-point2(2)).^2 + (point1(3)-point2(3)).^2);
%     l(2)=sqrt( (point1(1)-point3(1)).^2 + (point1(2)-point3(2)).^2 + (point1(3)-point3(3)).^2);
%     l(3)=sqrt( (point3(1)-point2(1)).^2 + (point3(2)-point2(2)).^2 + (point3(3)-point2(3)).^2);
%     s=sum(l)/2;
%     area(i)=real(sqrt(s*(s-l(1))*(s-l(2))*(s-l(3))));

end

end
