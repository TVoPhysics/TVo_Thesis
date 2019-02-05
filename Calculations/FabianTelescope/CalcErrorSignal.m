%%% Define the size and resolution
close all
clear all
dx=2e-5;
dy=dx;
x=[-1e-3:dx:1e-3]';
y=[-1e-3:dy:1e-3];

%%% Constants
z0=.1;
lambda=1.064e-6;

%%% Reference beam with a particular q-value
qx=1i*z0;
qy=qx;
xpos=1e-4;
vx0=.1e3;
inx=exp(-i*(pi/lambda/qx)*(x-xpos).^2)./sqrt(qx) .*exp(-i*2*pi*vx0*x);
iny=exp(-i*(pi/lambda/qy)*y.^2)./sqrt(qy);
inref=i*inx*iny;

mm = -.1:0.01:.1;

% Waist Location Shift
err_0deg_DeltaZ_list = [];
err_45deg_DeltaZ_list = [];

% Waist Size Shift
err_0deg_DeltaZr_list = [];
err_45deg_DeltaZr_list = [];

for m=-.1:0.01:.1

%%% Input beam with a slight mismatch
dd=z0*m/5; % defines the mode mismatch in beam position
%dd=0;
qx=1i*z0+dd;
qy=qx;
inx=exp(-i*(pi/lambda/qx)*(x-xpos).^2)./sqrt(qx) .*exp(-i*2*pi*vx0*x);
iny=exp(-i*(pi/lambda/qy)*y.^2)./sqrt(qy);
in=i*inx*iny*exp(-i*1/280.0001666665*dd/1e-4);

%%% Fourier transform the input beam only in the x direction
[out,dvx]=FourierX(in,dx);
hlf=(size(x,1)-1)/2;

%%%Inverse Fourier Transform
vx=dvx*[-hlf:1:hlf]';
[in2,dx2]=iFourierX(out,dvx);
fx=z0/(1+1/sqrt(2));
DX=sqrt(2)*fx;

%%%Define the optical parameters of the telescope
hf=getLensX(fx,lambda,x,y);
HD=getSpaceX(DX,lambda,vx,y);
HmD2=getSpaceX(-DX/2,lambda,vx,y);

%%%Propogate the two beams through the optical system and sample at
%%%Gouy Phase 0 and 45
outFabian_45deg1= applyFabian_45degGouy(in,dx,hf,HD,HmD2);
outref_45deg1   = applyFabian_45degGouy(inref,dx,hf,HD,HmD2);

outFabian_0deg1= applyFabian(in,dx,hf,HD,HmD2);
outref_0deg1   = applyFabian(inref,dx,hf,HD,HmD2);

%%%%Calc Error signal Integrate the signal power and subtract
[err_0deg_pringle, err_0deg_diagonal] = CalcErr(outFabian_0deg1,outref_0deg1);
[err_45deg_pringle, err_45deg_diagonal] = CalcErr(outFabian_45deg1,outref_45deg1);
 
err_0deg_DeltaZ_list = [err_0deg_DeltaZ_list ; err_0deg_pringle];
err_45deg_DeltaZ_list = [err_45deg_DeltaZ_list ; err_45deg_pringle];

end
for m=-.1:0.01:.1

%%% Input beam with a slight waist size mismatch
dd=z0*m; % defines the mode mismatch
qx=1i*z0+1i*dd;
qy=qx;
inx=exp(-i*(pi/lambda/qx)*(x-xpos).^2)./sqrt(qx) .*exp(-i*2*pi*vx0*x);
iny=exp(-i*(pi/lambda/qy)*y.^2)./sqrt(qy);
in=i*inx*iny*exp(-i*1/2000.0001666665*dd/1e-4);

%%% Fourier transform the input beam only in the x direction
[out,dvx]=FourierX(in,dx);
hlf=(size(x,1)-1)/2;

%%%Inverse Fourier Transform
vx=dvx*[-hlf:1:hlf]';
[in2,dx2]=iFourierX(out,dvx);
fx=z0/(1+1/sqrt(2));
DX=sqrt(2)*fx;

%%%Define the optical parameters of the telescope
hf=getLensX(fx,lambda,x,y);
HD=getSpaceX(DX,lambda,vx,y);
HmD2=getSpaceX(-DX/2,lambda,vx,y);

%%%Propogate the two beams through the optical system

outFabian_45deg= applyFabian_45degGouy(in,dx,hf,HD,HmD2);
outref_45deg   = applyFabian_45degGouy(inref,dx,hf,HD,HmD2);

outFabian_0deg= applyFabian(in,dx,hf,HD,HmD2);
outref_0deg   = applyFabian(inref,dx,hf,HD,HmD2);

%%%%Calc Error signal
[err_0deg_pringle, err_0deg_diagonal] = CalcErr(outFabian_0deg,outref_0deg);
[err_45deg_pringle, err_45deg_diagonal] = CalcErr(outFabian_45deg,outref_45deg);

err_0deg_DeltaZr_list = [err_0deg_DeltaZr_list ; err_0deg_pringle];
err_45deg_DeltaZr_list = [err_45deg_DeltaZr_list ; err_45deg_pringle];

end


figure(1)
p=plot(mm.',err_0deg_DeltaZ_list, mm.',err_45deg_DeltaZ_list);
p(1).LineWidth = 5;
p(2).LineWidth = 5;
lgd=legend('QPD+ModeConverter@0Deg','QPD+ModeConverter@45Deg')
lgd.FontSize = 14;
title('Shifted Beam Position')
xlabel('$\Delta z (m)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
grid on

figure(2)
p=plot(mm.',err_0deg_DeltaZr_list, mm.',err_45deg_DeltaZr_list);
p(1).LineWidth = 5;
p(2).LineWidth = 5;
lgd=legend('QPD+ModeConverter@0Deg','QPD+ModeConverter@45Deg');
lgd.FontSize = 14;
title('Shifted Beam Size')
xlabel('$\Delta z_R (m)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
grid on


%%%Load the BPD signals from python, BPD responses are in mW, must convert
%%%to counts.
load('./finesse_BPD_shifted_z.mat');
load('./finesse_BPD_shifted_zR.mat');

figure(3)
p=plot(Shifted_z, -real(BPD1_45deg)*64000, Shifted_z, -real(BPD2_45deg)*64000 );
p(1).LineWidth = 5;
p(2).LineWidth = 5;
lgd=legend('BPD1 @ 0-deg gouy ghase','BPD2 @ 45-deg gouy phase');
lgd.FontSize = 14;
xlabel('$\Delta z (m)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
xlim([-.1, .1])
grid on

figure(4)
p=plot(Shifted_zR, -real(BPD1_0deg)*64000, Shifted_zR, -real(BPD2_0deg)*64000 );
p(1).LineWidth = 5;
p(2).LineWidth = 5;
lgd=legend('BPD1 @ 0-deg gouy ghase','BPD2 @ 45-deg gouy phase');
lgd.FontSize = 14;
xlabel('$\Delta z_R (m)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
xlim([-.1, .1])
grid on

