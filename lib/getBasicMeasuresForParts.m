function [AvgThickness,TotalArea,SmoothArea,lblNames]=getBasicMeasuresForParts(label,labels_adjParts,pialf,pialv,opialf,opialv,Thickness,whitef,whitev)

%get the basic measures for a list of lables (as many entries as unique
%label names exist)

    lblNames=unique(label);
    AvgThickness=zeros(length(lblNames),1);
    TotalArea=zeros(length(lblNames),1);
    SmoothArea=zeros(length(lblNames),1);
    
    %make facebased thickness
    ThicknessFB=makeFacebased(Thickness,pialf);
    

    for k=1:length(lblNames)
        ids=find(label==lblNames(k));
        

        %find faceids
        [TotalAreai,fid]=calcPartAreai(pialf,pialv,ids);
        TotalArea(k)=sum(TotalAreai);

        %find the smooth area
        ids=find(labels_adjParts==lblNames(k));
        SmoothArea(k)=calcPartArea(opialf,opialv,ids);
        
        %avg thickness weighted by areas:
        pw=TotalAreai(fid>0)/TotalArea(k);
        [wTotalAreai,wfid]=calcPartAreai(whitef,whitev,ids);
        swTotalAreai=sum(wTotalAreai);
        ww=wTotalAreai(wfid>0)/swTotalAreai;
        

        
        AvgThickness(k)=(sum(ThicknessFB(fid>0).*pw)+sum(ThicknessFB(wfid>0).*ww))/2;
       
    end