function [Wi,W]=deriveW(eta,psi,mdim,akz)
%%
mdim2=mdim;
[NROW,NCOL]=size(eta);
z0=0;
phi=zeros(NROW,NCOL,mdim2);
phizz0=psi;
phi(:,:,1)=psi;



%ê€ìÆìWäJ
for N=2:mdim2

    for k=1:N-1
        phi(:,:,N)=phi(:,:,N)-...
           (eta-z0).^k/factorial(k).*real(ifft(akz.^k.*fft(phi(:,:,N-k),NCOL,2),NCOL,2));

    %(eta-z0).^k/factorial(k).*real(fft1d2d(akz.^k.*fft1d2d(phi(:,:,N-k),2,1),2,-1));

    end
    phizz0=phizz0+phi(:,:,N);

end


W=zeros(NROW,NCOL);
Wi=zeros(NROW,NCOL,mdim2);

for N=1:mdim2
    for k=1:N;
        Wi(:,:,N)=Wi(:,:,N)+(eta).^(k-1)/factorial(k-1).*real(ifft(akz.^k.*fft(phi(:,:,N-k+1),NCOL,2),NCOL,2));
%        Wi(:,:,N)=Wi(:,:,N)+(eta).^(k-1)/factorial(k-1).*real(fft1d2d(akz.^k.*fft1d2d(phi(:,:,N-k+1),2,1),2,-1));
    end
    W=W+Wi(:,:,N);

end
    

