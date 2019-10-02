function [d,varargout]=cohensD(data,grp)
    cats=unique(grp);
    
    data1=data(grp==cats(1));
    data2=data(grp==cats(2));
    
    n1=numel(data1);
    n2=numel(data2);
    
    m1=mean(data1);
    m2=mean(data2);
    
    s1=std(data1);
    s2=std(data2);
    
    sp=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
    d=(m2-m1)/sp;
    
    
    if nargout==2%calculate p value?
        p = ranksum(data1,data2);
        varargout{1} = p;
    end

end