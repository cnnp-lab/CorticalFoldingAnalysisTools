function [totalPialArea,totalPialSmoothArea,totalAvgThickness]=extract_PialArea_SmoothArea_AvgThickness(Thickness,pialv,pialf,opialv,opialf)
    Thicknessf=makeFacebased(Thickness,pialf);
    
    pialarea=calcTriangleArea(pialf,pialv);
    opialarea=calcTriangleArea(opialf,opialv);
    
    faceids=find(Thicknessf==0);
    
    idssmooth=findcorrsmooth(faceids,opialf,opialv,pialf,pialv);
    %c=zeros(length(opialf),1); c(idssmooth)=1;
    %trisurf(opialf,opialv(:,1),opialv(:,2),opialv(:,3),c)
    
    foldedIDs=1:length(pialf);
    foldedIDs(faceids)=[];
    smoothIDs=1:length(opialf);
    smoothIDs(idssmooth)=[];
    
    totalPialArea=sum(pialarea(foldedIDs));
    totalPialSmoothArea=sum(opialarea(smoothIDs));
    totalAvgThickness=mean(Thicknessf(foldedIDs));