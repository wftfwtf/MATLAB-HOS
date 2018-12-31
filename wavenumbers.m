function [varargout]=wavenumbers(xdomain,ydomain,NCOL,NROW)
 dkx=2*pi/xdomain;
 dky=2*pi/ydomain;
 akx=zeros(1,NCOL);
 aky=zeros(1,NROW);
akz=zeros(NROW,NCOL);
 
for i=1:NCOL
    if i <= ceil(NCOL/2)
        akx(i)=dkx*(i-1);
    else
        akx(i)=dkx*(i-NCOL-1);
    end
end

for i=1:NROW
    if i<= ceil(NROW/2)
        aky(i)=dky*(i-1);
    else
        aky(i)=dky*(i-NROW-1);
    end
end

for j=1:NCOL
    for i=1:NROW
        akz(i,j)=sqrt(akx(j)^2+aky(i)^2);
    end
end

varargout{1}=akx;
varargout{2}=aky;
varargout{3}=akz;

