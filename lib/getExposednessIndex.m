function [EI,nBB,EIf,varargout]=getExposednessIndex(in_exposed,varargin)
    SE = strel('sphere',1);
    K=SE.Neighborhood;
    K(2,2,2)=0;
    C = convn(in_exposed,K,'same');
    mSE=sum(K(:));
    C=mSE-C;
    C(C==mSE)=0;
    C(in_exposed==0)=0;
%     imagesc(squeeze(C(100,:,:)))
    EI=sum(C(:));%
    Ci=C>0;
    nBB=sum(Ci(:));%number of boundary boxes
    
    C(C==1)=4/3;
    C(C==2)=4/3;
    C(C==3)=1;
    C(C==4)=2/3;
    C(C==5)=sqrt(5)/3;
    EIf=sum(C(:));
    
    if nargout>3 %calculate more stuff from voxelisation
        scale=varargin{1};
        %calculate slice perimeters in x-y direction
        K=zeros(3,3,1);
        K(:,:,1)=[0 1 0;
                  1 0 1;
                  0 1 0];
        s=sum(K(:));
        C = convn(in_exposed,K,'same');

        C(in_exposed==0)=NaN;
        C(C==s)=NaN;

        C(C==1)=(1+sqrt(2))*scale;
        C(C==2)=(sqrt(2))*scale;
        C(C==3)=scale;

        %perimeter: 
        stackedperimeter=squeeze(nansum(nansum(C,1),2));
        %area:
        stackedarea=squeeze(sum(sum(in_exposed,1),2).*scale^2);
        
        %estimate surface area:
        estimatedarea=sum(sqrt((stackedarea(1:end-1)-stackedarea(2:end)).^2+scale^2.*(stackedperimeter(1:end-1)+stackedperimeter(2:end)).^2/4));
        %estimate volume:
        estimatedvol=sum(scale/3*(stackedarea(1:end-1)+stackedarea(2:end)+sqrt(stackedarea(1:end-1).*stackedarea(2:end))));
        
        varargout{1}=estimatedarea;
        varargout{2}=estimatedvol;
    end
    
    
end