function olabel=matchSurfLabel(label,pialv,opialv)

olabel=zeros(length(opialv),1);
for k=1:length(opialv)
    cp=opialv(k,:);
%     d=sqrt((pialv(:,1)-cp(1)).^2+(pialv(:,2)-cp(2)).^2+(pialv(:,3)-cp(3)).^2);
    d=sqrt(sum((pialv-cp).^2,2));
    [~,mid]=min(d);
    olabel(k)=label(mid);
end