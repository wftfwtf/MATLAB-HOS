function [varargout]=compHOSMRK4_MATLAB(eta0,psi0,mdim,xdomain,ydomain,...
    grav,deltat,NT,Ta,istep_out,Nx,xp,dispflag,varargin)


restart=1;
if nargin>13
    restart=varargin{1};
end

[r,c]=size(eta0);

if r==Nx;
    eta0=eta0';psi0=psi0';
    Nens=c;
else
    Nens=r;
end

obsrec=zeros(Nens,round((NT-restart)/istep_out));
ETAPSIa=zeros(round((NT-restart)/istep_out),2*Nx);

conservQuan=zeros(round((NT-restart)/istep_out),3);

dx=xdomain/Nx;

Hs=4*std(eta0(1,:));
Hspsi=4*std(psi0(1,:));

if sign(deltat)>0
    istepfirst=restart+1;isteplast=NT;
else
    istepfirst=NT;isteplast=restart+1;
end

for istep=istepfirst:sign(deltat):isteplast
    % W is for buoy time records
    [eta0,psi0,W]=compHOSMRK4(eta0,psi0,mdim,xdomain,ydomain,grav,deltat,istep,Ta,1);
    
    if mod(istep,istep_out)==0
        
        ti=round((istep-restart)/istep_out);
        obsrec(:,ti)=eta0(:,round(xp/dx));
        
        if nargout>1
            ETAPSIa(ti,:)=[eta0(1,:) psi0(1,:)];
        end
        
        if nargout>2
            conservQuan(ti,:)=[sum(eta0(1,:)) sum(eta0(1,:).^2) sum(psi0(1,:).*W(1,:))];
        end
        
        if  dispflag
            istep
            if dispflag>1
                figure(11);
                subplot(2,1,1);
                plot(eta0');title(num2str(istep));ylim([-2*Hs 2*Hs])
                subplot(2,1,2);
                plot(psi0');title(num2str(istep));ylim([-2*Hspsi 2*Hspsi])
                drawnow
            end
        end
    end
end

varargout{1}=obsrec;

if nargout>1
    varargout{2}=ETAPSIa;
    if nargout>2
        varargout{3}=conservQuan;
    end
end

end