function [weightedAvgThicknessGM,weightedAvgThicknessWM]=extract_weightedAvgThickness(Thickness,pialv,pialf,whitev,whitef)
    Thicknessf=makeFacebased(Thickness,pialf);
    
    pialarea=calcTriangleArea(pialf,pialv);
    npialarea=pialarea/sum(pialarea);
    whitearea=calcTriangleArea(whitef,whitev);
    nwhitearea=whitearea/sum(whitearea);
    %opialarea=calcTriangleArea(opialf,opialv);
    
    faceids=find(Thicknessf==0);
    
    
    foldedIDs=1:length(pialf);
    foldedIDs(faceids)=[];

    weightedAvgThicknessGM=sum(Thicknessf(foldedIDs).*npialarea(foldedIDs));
    weightedAvgThicknessWM=sum(Thicknessf(foldedIDs).*nwhitearea(foldedIDs));
end