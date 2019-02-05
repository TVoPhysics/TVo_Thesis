function [ ] = movieHGconverter( m,n,step )
%Plots and saves images procued by passing a mnHG mode through a the
% mode converting telescope and rotates the input from 0 to 360 degrees in 

%% Example of usage 
% movieHGcoverter(1,1,10)
% produces a 11HG mode and rotates it from 0 to 360 in 10 steps.
mkdir 'movieHGconverterImages'
cd ../..
addpath software/matlab/HG/
addpath AMM_mode_converter_QPD_sensor/FabianTelescope/

for count=0:step
    %% Define some constant parameters and point to location of functions
    addpath 'Mode_converter_QPD/Fabian_study/HG/'
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
    %m = 1;                                                                      % Choose m value of HG mode
    %n = 1;                                                                      % Choose n value of HG mode
    angle=count*floor(360/step);                                                 % Choose telescope rotation angle
    inFabian = rotatedHGmode(m,n,x,y,zx,zRx,zy,zRy,lambda,angle);

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

    outFabian=applyFabian(inFabian,dx,hf,HD,HmD2);
    
    % Plot the gaussian beam with mode converter graphic
    plotModeConverter(inFabian,outFabian);
    orient landscape;
    saveas(gcf,['AMM_mode_converter_QPD_sensor/FabianTelescope/movieHGconverterImages/angle',num2str(angle),'.png']);
    close all
end

cd AMM_mode_converter_QPD_sensor/FabianTelescope/

end

