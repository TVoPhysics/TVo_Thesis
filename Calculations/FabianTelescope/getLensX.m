function [h]=getLensX(fx,lambda,xaxis,yaxis)

xx=xaxis*ones(size(yaxis));

h=exp(i*(pi/lambda/fx)*xx.^2);

%contour(xaxis*ones(size(yaxis)),ones(size(xaxis))*yaxis,real(h))
%colorbar
