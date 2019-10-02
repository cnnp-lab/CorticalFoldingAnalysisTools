function [totalPialArea,totalWhiteArea,totalPialSmoothArea,totalAvgThickness]=extract_PialArea_WhiteArea_SmoothArea_AvgThickness(Thickness,pialv,pialf,whitev,whitef,opialv,opialf)
    Thicknessf=makeFacebased(Thickness,pialf);
    
    pialarea=calcTriangleArea(pialf,pialv);
    whitearea=calcTriangleArea(whitef,whitev);
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
    totalWhiteArea=sum(whitearea(foldedIDs));
    totalPialSmoothArea=sum(opialarea(smoothIDs));
    totalAvgThickness=mean(Thicknessf(foldedIDs));