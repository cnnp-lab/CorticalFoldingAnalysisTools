function [mslope,stdslope,mpslope,moffset,stdoffset,mpoffset] = estimateRegressionConfidence_twogroupdiff(xs,ys,leaven,mtimes)

slope=zeros(mtimes,2);
offset=slope;
p_slope=slope;
p_offset=slope;

for m=1:mtimes
    ids=randi(length(xs),leaven,1);
    
    xc=xs;
    yc=ys;
    
    xg1=xc(ids);
    yg1=yc(ids);
    
    xc(ids)=[];
    yc(ids)=[];
    
    crt=table(xg1,yg1);
    fitsg=fitlm(crt,'yg1~xg1');

    slope(m,1)=table2array(fitsg.Coefficients(2,1));
    offset(m,1)=table2array(fitsg.Coefficients(1,1));

    p_slope(m,1)=table2array(fitsg.Coefficients(2,4));
    p_offset(m,1)=table2array(fitsg.Coefficients(1,4));
    
    
    
    
    ids=randi(length(xc),leaven,1);
    
    xg2=xc(ids);
    yg2=yc(ids);
    
    
    crt=table(xg2,yg2);
    fitsg=fitlm(crt,'yg2~xg2');

    slope(m,2)=table2array(fitsg.Coefficients(2,1));
    offset(m,2)=table2array(fitsg.Coefficients(1,1));

    p_slope(m,2)=table2array(fitsg.Coefficients(2,4));
    p_offset(m,2)=table2array(fitsg.Coefficients(1,4));
    

end

slope(p_slope>0.01)=[];
offset(p_slope>0.01)=[];
p_offset(p_slope>0.01)=[];
p_slope(p_slope>0.01)=[];

mslope=median(slope(:));
stdslope=std(slope(:,1)-slope(:,2));
mpslope=median(p_slope(:));

moffset=median(offset(:));
stdoffset=std(offset(:,1)-offset(:,2));
mpoffset=median(p_offset(:));


