%% Setting


kp=1;grav=1;omegap=sqrt(grav*kp); % normalize so that grav=1
steepness=0.11;
Hm0=2*steepness/kp;gamma=3.3;
mdim=3;
xdomain=128*2*pi/kp;ydomain=1; NCOL=(xdomain/2/pi*kp)*8*(mdim+1);NROW=1;

xp=xdomain/2;
[akx,aky,AKZ]=wavenumbers(xdomain,ydomain,NCOL,NROW);
afilter=antialias(NROW,NCOL,mdim,0);
omega=sqrt(grav*AKZ);
omega(1)=Inf;

dx=xdomain/NCOL;
dkx=2*pi/xdomain;
X=(1:NCOL)*dx;

dt_per_Tm=80; % dt/Tm=80
Nperiod=50; % 50 wave periods
Taperiod=5;
Ta=Taperiod/2.5*dt_per_Tm;

deltat=(2*pi/omegap)/dt_per_Tm;

NT=(Nperiod+Taperiod)*dt_per_Tm;


w=linspace(0,max(akx),10000);

S=jonswap(w,[Hm0 2*pi/omegap 3.3 0.07 0.09 -1],0);
S.g=grav;
Sk=time2spa(S,akx(1:end/2),[],grav);

Sk2=Sk;
Sk2.S=zeros(1,NCOL);Sk2.S(1:length(Sk.S))=Sk.S;
Sk2.k=zeros(1,NCOL);Sk2.k=akx;

%% Time evolution and detect freak wave

maxEta=0

%while maxEta<1.3*Hm0

feta0=sqrt(2*Sk2.S*dkx).*exp(1i*2*pi*rand(size(Sk2.S)))*NCOL.*afilter;
eta0=real(ifft(feta0));
psi0=real(ifft(1i*grav./omega.*feta0));

(4*std(eta0)-Hm0)/Hm0

figure(4)
subplot(2,1,1)
plot(X/2/pi,eta0)
subplot(2,1,2)
plot(X/2/pi,psi0)
drawnow

%%
Tout_per_Tm=20;
istep_out=dt_per_Tm/Tout_per_Tm;
tic;
[obsrec,ETAPSIt,conservQuan]=compHOSMRK4_MATLAB(eta0,psi0,mdim,xdomain,ydomain,...
    grav,deltat,NT,Ta,istep_out,NCOL,xp,2);
toc

maxEta=max(max(ETAPSIt(:,1:end/2)));


maxEta/Hm0
%end

%% Plot

figure(6)
pcolor(ETAPSIt(:,1:end/2))
shading flat

figure(7)
tmp=sum(conservQuan(:,2:3),2);
plot((tmp/tmp(1)-1)*100)
