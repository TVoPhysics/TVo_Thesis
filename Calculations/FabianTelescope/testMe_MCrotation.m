%% Define some constant parameters and point to location of functions
cd ../..
addpath software/matlab/HG/
addpath AMM_mode_converter_QPD_sensor/FabianTelescope/
lambda=1.064e-6;
dx=2e-6;
dy=dx;

x=[-1e-3:dx:1e-3]';
y=[-1e-3:dy:1e-3];

z0=1e-1;
qx=1i*z0;
qy=qx;

zx=0;
zy=zx;
zRx=imag(qx);
zRy=zRx;

%% Using getHG calculate inFabian
m = 1;                                                                      % Choose m value of HG mode
n = 1;                                                                      % Choose n value of HG mode
angle=-22.5;                                                                   % Choose telescope rotation angle
%unnormal = rotatedHGmode(m,n,x,y,zx,zRx,zy,zRy,lambda,angle);
%unmax = max(max(unnormal));
%unmin = min(min(unnormal));
inFabian = rotatedHGmode(m,n,x,y,zx,zRx,zy,zRy,lambda,angle);%(unnormal - unmin)./(unmax - unmin);

%% Applying the mode converter telescope with Fourier optics
[out,dvx]=FourierX(inFabian,dx);
hlf=(size(x,1)-1)/2;
vx=dvx*[-hlf:1:hlf]';
[inFabian2,dx2]=iFourierX(out,dvx);

fx=z0/(1+1/sqrt(2));
DX=sqrt(2)*fx;
hf=getLensX(fx,lambda,x,y);
HD=getSpaceX(DX,lambda,vx,y);
HmD2=getSpaceX(-DX/2,lambda,vx,y);

applyFabian(inFabian,dx,hf,HD,HmD2);
outFabian = applyFabian(inFabian,dx,hf,HD,HmD2);

% Plot the gaussian beam with mode converter graphic
plotModeConverter(inFabian,outFabian);
cd AMM_mode_converter_QPD_sensor/FabianTelescope/
