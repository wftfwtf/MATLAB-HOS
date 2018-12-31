# MATLAB-HOS
MATLAB version of HOS (Higher Order Spectral method; Dommermuth & Yue 1987, West et al. 1987).

        [obsrec,ETAPSIa,conservQuan]=compHOSMRK4_MATLAB(etatmp,psitmp,mdim,xdomain,ydomain,...
 
        grav,deltat,NT,Ta,istep_out,NCOL,xp,dispflag);

obsrec: Observation record at a point

ETAPSIa: Simulated surface elevation and velocity potential at the surface

conservQuan: Conserved Quantities

etatmp: Initial value of surface elevation

psitmp: Initial value of velocity potential at the surface

mdim: Nonlinearity order M

xdomain, ydomain:　Domain length

grav: Gravity accerelation

deltat: Interval of Time step

NT: Number of time steps

istep_out: Interval of output steps

NCOL: Number of spatial grids

xp: Cordinate of the observation point

dispflag: Animation displayed if >1
