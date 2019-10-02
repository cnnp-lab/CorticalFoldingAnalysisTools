function [slopeqtl,offsetqtl] = estimateRegressionConfidence(xs,ys,qtl,leaven,mtimes)

slope=zeros(mtimes,1);
offset=slope;
p_slope=slope;
p_offset=slope;

for m=1:mtimes
    ids=randi(length(xs),leaven,1);
    
    xc=xs;
    yc=ys;
    
    xc(ids)=[];
    yc(ids)=[];
    
    crt=table(xc,yc);
    fitsg=fitlm(crt,'yc~xc');

    slope(m)=table2array(fitsg.Coefficients(2,1));
    offset(m)=table2array(fitsg.Coefficients(1,1));

    p_slope(m)=table2array(fitsg.Coefficients(2,4));
    p_offset(m)=table2array(fitsg.Coefficients(1,4));

end

slope(p_slope>0.01)=[];
offset(p_slope>0.01)=[];
p_offset(p_slope>0.01)=[];
p_slope(p_slope>0.01)=[];


slopeqtl=quantile(slope,qtl);
offsetqtl=quantile(offset,qtl);
