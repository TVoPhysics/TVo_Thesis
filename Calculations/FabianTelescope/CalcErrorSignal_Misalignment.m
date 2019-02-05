%%% Define the size and resolution
dx=2e-5;
dy=dx;
x=[-1e-3:dx:1e-3]';
y=[-1e-3:dy:1e-3];

%%% Constants
z0=.1;
lambda=1.064e-6;

%%% Reference beam with a particular q-value
xpos=0e-4;
vx0=0;
w0 = sqrt(lambda*z0/pi);

inx=exp(-i*((x-xpos).^2)/w0^2).*exp(-i*2*pi*vx0*x);
iny=exp(-i*((y).^2)/w0^2);
inref=i*inx*iny;

mm = -.1:0.01:.1;

% Cavity Tilt Shift
err_0deg_DeltaA_list = [];
err_90deg_DeltaA_list = [];

% Cavity Pos Shift
err_0deg_DeltaX_list = [];
err_90deg_DeltaX_list = [];

for m=-.1:0.01:.1

%%% Input beam with a slight mismatch
% defines the mode mismatch in beam position
xpos1=1e-5*m;
vx0=0;
inx=exp(-i*((x-xpos).^2)/w0^2).*exp(-i*2*pi*vx0*x);
iny=exp(-i*((y).^2)/w0^2);
in=i*inx*iny;

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
outFabian_90deg1= applyFabian_90degGouy(in,dx,hf,HD,HmD2);
outref_90deg1   = applyFabian_90degGouy(inref,dx,hf,HD,HmD2);

outFabian_0deg1= applyFabian(in,dx,hf,HD,HmD2);
outref_0deg1   = applyFabian(inref,dx,hf,HD,HmD2);

%%%%Calc Error signal Integrate the signal power and subtract
[err_0deg_pringle, err_0deg_diagonal] = CalcErr(outFabian_0deg1,outref_0deg1);
[err_90deg_pringle, err_90deg_diagonal] = CalcErr(outFabian_90deg1,outref_90deg1);

err_0deg_DeltaX_list = [err_0deg_DeltaX_list ; err_0deg_pringle];
err_90deg_DeltaX_list = [err_90deg_DeltaX_list ; err_90deg_pringle];

end
for m=-.1:0.01:.1

%%% Input beam with a slight waist size mismatch
%dd=z0*m; % defines the mode mismatch
xpos1=0;
vx0=m;
inx=exp(-i*((x-xpos).^2)/w0^2).*exp(-i*2*pi*vx0*x);
iny=exp(-i*((y).^2)/w0^2);
in=i*inx*iny;

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

outFabian_90deg= applyFabian_90degGouy(in,dx,hf,HD,HmD2);
outref_90deg   = applyFabian_90degGouy(inref,dx,hf,HD,HmD2);

outFabian_0deg= applyFabian(in,dx,hf,HD,HmD2);
outref_0deg   = applyFabian(inref,dx,hf,HD,HmD2);

%%%%Calc Error signal
[err_0deg_pringle, err_0deg_diagonal] = CalcErr(outFabian_0deg,outref_0deg);
[err_90deg_pringle, err_90deg_diagonal] = CalcErr(outFabian_90deg,outref_90deg);

err_0deg_DeltaA_list = [err_0deg_DeltaA_list ; err_0deg_pringle];
err_90deg_DeltaA_list = [err_90deg_DeltaA_list ; err_90deg_pringle];

end

figure(1)
plot(mm.',err_0deg_DeltaX_list, mm.',err_90deg_DeltaX_list)
legend('QPD+ModeConverter@0Deg','QPD+ModeConverter@90Deg')
title('Shifted Beam Axis')
xlabel('$\Delta x_{pos} (m)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
grid on

figure(2)
plot(mm.',err_0deg_DeltaA_list, mm.',err_90deg_DeltaA_list)
legend('QPD+ModeConverter@0Deg','QPD+ModeConverter@90Deg')
title('Shifted Beam Angle')
xlabel('$\Delta \alpha_x (deg)$','interpreter','latex','FontSize',18')
ylabel('Response (Counts)','interpreter','latex','FontSize',18')
grid on

figure(3)
contour(abs(outFabian_0deg))