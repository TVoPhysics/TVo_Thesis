function [out,dvx]=FourierX(in,dx)
% x is first dimension
% out = int in(x) * exp(i*2*pi*nux*x) dx

d=1;
N=size(in,d);
dvx=1/N/dx;

out=N*dx*fftshift(ifft(ifftshift(in,d),N,d),d);

