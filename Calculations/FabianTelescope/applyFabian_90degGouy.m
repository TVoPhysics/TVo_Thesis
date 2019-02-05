function [out]=applyFabian_90degGouy(in,dx,hf,HD,HmD2)

[FD,dvx]=FourierX(in,dx);
FD=HmD2.*FD;
XD=iFourierX(FD,dvx);
XD=hf.*XD;
FD=FourierX(XD,dx);
FD=HD.*FD;
XD=iFourierX(FD,dvx);
XD=hf.*XD;
FD=FourierX(XD,dx);
FD=10*HD.*FD;
XD=iFourierX(FD,dvx);
out=XD/sqrt(i);