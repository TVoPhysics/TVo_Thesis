%Inner radius
R0=1e-3/2;

%Outer radius
Rm=5e-3/2;

% Define BPD array
dx=0.00001;
dy=dx;
edge=0.0025;
x=(-edge:dx:edge);
y=(-edge:dy:edge)';
X=ones(size(y))*x;
Y=y*ones(size(x));


M1=X.^2+Y.^2<R0^2;
Ri=and(X.^2+Y.^2>=R0^2,X.^2+Y.^2<Rm^2);
Sp=180/pi*angle(X+i*Y);
Sq=mod(180/pi*angle(X+i*Y),360);

W2=and(Sp>=330,Sp<90);
W3=and(Sq>=90,Sq<210);
W4=and(Sp>=210,Sp<330);

M2=and(W2,Ri);
M3=and(W3,Ri);
M4=and(W4,Ri);

csvwrite('center.dat',M1)
csvwrite('upperright.dat',M2)
csvwrite('upperleft.dat',M3)
csvwrite('lowercenter.dat',M2)