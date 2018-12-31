function [dteta,dtpsi,varargout]=compHOSM(eta0,psi0,mdim,xdomain,ydomain,grav,istep,Ta,dim)

[NROWens,NCOL]=size(eta0);

if dim==1
    NROW=1;ens=NROWens;
else %dim==2
    NROW=NROWens;ens=1;
end

[akx,aky,AKZ]=wavenumbers(xdomain,ydomain,NCOL,NROW);
[AKX,AKY]=meshgrid(akx,aky);
afilter=antialias(NROW,NCOL,mdim,0);

if dim==1
    AKX=repmat(AKX,ens,1);AKZ=repmat(AKZ,ens,1);
    afilter=repmat(afilter,ens,1);
end
[Wi,W]=deriveW(eta0,psi0,mdim,AKZ);

dteta=zeros(NROWens,NCOL);dtpsi=zeros(NROWens,NCOL);
dyeta=zeros(NROWens,NCOL);dypsi=zeros(NROWens,NCOL);

if mdim>=2 % Nonlinear
        dxeta=real(ifft(1i*AKX.*fft(eta0,NCOL,2),NCOL,2));
        dxpsi=real(ifft(1i*AKX.*fft(psi0,NCOL,2),NCOL,2));

    
    if dim==2
        dyeta=real(fft1d2d(1i*AKY.*fft1d2d(eta0,dim,1),dim,-1));
        dypsi=real(fft1d2d(1i*AKY.*fft1d2d(psi0,dim,1),dim,-1));
    end
    
    greta2=dxeta.^2+dyeta.^2;
    
    dteta=-dxpsi.*dxeta-dypsi.*dyeta+Wi(:,:,2);
    dtpsi=-1/2*(dxpsi.^2+dypsi.^2)+1/2*Wi(:,:,1).^2;
    
    if mdim>=3 % over than 3rd order
        for imdim=3:mdim;
            for i=3:imdim
                dteta=dteta+greta2.*Wi(:,:,i-2)+Wi(:,:,i);
            end
            for i=1:imdim-1;
                dtpsi=dtpsi+1/2*(Wi(:,:,i).*Wi(:,:,imdim-i));
            end
            
            if imdim>=4
                for i=1:imdim-3;
                    dtpsi=dtpsi+1/2*(Wi(:,:,i).*Wi(:,:,imdim-2-i)).*greta2;
                end
            end
        end
    end
end

F=1-exp(-(istep/Ta)^2);
dtpsi=dtpsi*F-grav*eta0;
dteta=dteta*F+Wi(:,:,1);

% Anti-aliasing

dtpsi=real(ifft(afilter.*fft(dtpsi,NCOL,2),NCOL,2));
dteta=real(ifft(afilter.*fft(dteta,NCOL,2),NCOL,2));
% dtpsi=real(fft1d2d(afilter.*fft1d2d(dtpsi,dim,1),dim,-1));
% dteta=real(fft1d2d(afilter.*fft1d2d(dteta,dim,1),dim,-1));

if nargout>2
    varargout{1}=W;
end