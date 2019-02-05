function [H]=getSpaceX(DX,lambda,vxaxis,yaxis)

vxx=vxaxis*ones(size(yaxis));

H=exp(i*(pi*lambda*DX)*vxx.^2);

%contour(xaxis*ones(size(yaxis)),ones(size(xaxis))*yaxis,real(h))
%colorbar
