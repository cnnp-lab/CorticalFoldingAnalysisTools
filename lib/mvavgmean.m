function [tmv,ampl,stdv]=mvavgmean(t,x,wl)
    tstart=min(t); tend=max(t);
    nwindow=floor((tend-tstart)/wl);
    ampl=zeros(1,nwindow);
    tmv=zeros(1,nwindow);
    stdv=zeros(1,nwindow);
    for k=1:nwindow
        ts=(k-1)*wl+tstart;
        te=k*wl+tstart;
        xs=x(t>=ts&t<te);
        tmv(k)=(k-0.5)*wl+tstart;
        ampl(k)=mean(xs);
        stdv(k)=std(xs);
    end
end