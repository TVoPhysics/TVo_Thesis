dx=2e-6;
dy=dx;
x=[-1e-3:dx:1e-3]';
y=[-1e-3:dy:1e-3];

z0=1e-1;
lambda=1.064e-6;

qx=1i*z0;
qy=qx;
xpos=0e-4;
vx0=0e3;
inx=exp(-i*(pi/lambda/qx)*(x-xpos).^2)./sqrt(qx) .*exp(-i*2*pi*vx0*x);
iny=exp(-i*(pi/lambda/qy)*y.^2)./sqrt(qy);
inref=i*inx*iny;

dd=z0*0.001; % defines the mode mismatch
qx=1i*z0+dd;
qy=qx;
inx=exp(-i*(pi/lambda/qx)*(x-xpos).^2)./sqrt(qx) .*exp(-i*2*pi*vx0*x);
iny=exp(-i*(pi/lambda/qy)*y.^2)./sqrt(qy);
in=i*inx*iny*exp(-i*1/2000.0001666665*dd/1e-4);


[out,dvx]=FourierX(in,dx);
hlf=(size(x,1)-1)/2;
vx=dvx*[-hlf:1:hlf]';
[in2,dx2]=iFourierX(out,dvx);

fx=z0/(1+1/sqrt(2));
DX=sqrt(2)*fx;
hf=getLensX(fx,lambda,x,y);
HD=getSpaceX(DX,lambda,vx,y);
HmD2=getSpaceX(-DX/2,lambda,vx,y);

outFabian=applyFabian(in,dx,hf,HD,HmD2);
outref   =applyFabian(inref,dx,hf,HD,HmD2);

%plot(x,real(in),x,real(in2));
%contour(x*ones(size(y)),ones(size(x))*y,(abs(in)));
%contour(vx*ones(size(y)),ones(size(x))*y,(abs(out)));

figure(10)
contour(x*ones(size(y)),ones(size(x))*y,(abs(inref)));
figure(11)
contour(x*ones(size(y)),ones(size(x))*y,(abs(outref)));
figure(2)
contour(x*ones(size(y)),ones(size(x))*y,(imag( outFabian.*conj(outref) )));
colorbar
figure(3)
contour(x*ones(size(y)),ones(size(x))*y,(imag( in.*conj(inref) ))); % *exp(-i/2000)
colorbar
sum(sum(imag( in.*conj(inref))))
sum(sum(imag( outFabian.*conj(outref))))
