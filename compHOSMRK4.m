function [eta2,psi2,W]=compHOSMRK4(eta0,psi0,mdim,xdomain,ydomain,grav,deltat,istep,Ta,dim)

    eta=eta0;psi=psi0;


    [dteta,dtpsi]=compHOSM(eta,psi,mdim,xdomain,ydomain,grav,istep,Ta,dim);
    rk1e=dteta;rk1p=dtpsi;
    eta=eta0+dteta*deltat/2;psi=psi0+dtpsi*deltat/2;
    
    % 2nd RK4
    [dteta,dtpsi]=compHOSM(eta,psi,mdim,xdomain,ydomain,grav,istep,Ta,dim);
    rk2e=dteta;rk2p=dtpsi;
    eta=eta0+dteta*deltat/2;psi=psi0+dtpsi*deltat/2;
    
    % 3rd RK4
    [dteta,dtpsi]=compHOSM(eta,psi,mdim,xdomain,ydomain,grav,istep,Ta,dim);
    rk3e=dteta;rk3p=dtpsi;
    eta=eta0+dteta*deltat;psi=psi0+dtpsi*deltat;
    
    % 4th RK4
    [dteta,dtpsi,W]=compHOSM(eta,psi,mdim,xdomain,ydomain,grav,istep,Ta,dim);
    rk4e=dteta;rk4p=dtpsi;
    eta2=eta0+deltat/6*(rk1e+2*rk2e+2*rk3e+rk4e);
    psi2=psi0+deltat/6*(rk1p+2*rk2p+2*rk3p+rk4p);