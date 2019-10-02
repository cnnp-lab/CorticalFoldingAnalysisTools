function d=cohensD2inp(data1,data2)
    
    
    n1=numel(data1);
    n2=numel(data2);
    
    m1=mean(data1);
    m2=mean(data2);
    
    s1=std(data1);
    s2=std(data2);
    
    sp=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
    d=(m2-m1)/sp;

end