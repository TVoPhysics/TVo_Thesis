function [out,dx]=iFourierX(in,dvx)
% x is first dimension
% out = int in(x) * exp(-i*2*pi*nux*x) dx

d=1;
N=size(in,d);
dx=1/N/dvx;

out=dvx*fftshift(fft(ifftshift(in,d),N,d),d);