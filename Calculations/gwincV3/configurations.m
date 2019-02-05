function [y,f,sc] = configurations()
%startup;

%Optimized for NSNS:
%[input power, SRC phase, T_SRM, T_ITM, T_PRM, homodyne phase]

NOSRM = [25,   0,          1, 0.014, 0.027, 130/180*pi]; 
NSNS=   [125, 11/180*pi, 0.2, 0.014, 0.027, 103/180*pi];
ZERODET=[125,  0       , 0.2, 0.014, 0.027, 116/180*pi];
ZDETLP =[ 25,  0       , 0.2, 0.014, 0.027, 116/180*pi];
NSNOSRM=[ 40,  0,          1, 0.014, 0.027, 120/180*pi];

% Optimized for BHBH
BHBH=   [4.5, 72/180*pi, 0.2, 0.014, 0.027,  90/180*pi];
BHBH20= [ 20, 20/180*pi, 0.2, 0.014, 0.027, 105/180*pi];
BHBH30= [ 19, 30/180*pi, 0.2, 0.014, 0.027,  78/180*pi];
BHNOSRM=[7.8,  0       ,   1, 0.014, 0.027, 145/180*pi];

% optimized for 1kHz
KHZ=    [125,4.7/180*pi,0.011,0.014, 0.027, 128/180*pi];

EXTREMEHF=[125,63.4/180*pi,0.03,0.21,0.002,135/180*pi];


%conjGradientOptimizer(0,NSNS,[0,0,0,0,1,0])
[y(1,:),f,sc(1,:)]=callGwinc(NOSRM);
[y(2,:),f,sc(2,:)]=callGwinc(ZDETLP);
[y(3,:),f,sc(3,:)]=callGwinc(ZERODET);
[y(4,:),f,sc(4,:)]=callGwinc(NSNS);
[y(5,:),f,sc(5,:)]=callGwinc(BHBH20);
[y(6,:),f,sc(6,:)]=callGwinc(BHBH);
[y(7,:),f,sc(7,:)]=callGwinc(KHZ);

set(0,'DefaultAxesLineStyleOrder',{'-','-.*','-o'})
%hndls=loglog(f,sqrt(y));
hndls=loglog(f,sqrt(y((1:5),:)));

grid on;
%xlabel('Hz');
%ylabel('strain 1/rtHz');
%title('AdvLIGO tunings');
%legend('0) NO SRM','1) Zero Detune, low power','2) Zero Detune, high power','NS-NS Opt.',...
%      'BH-BH 20 deg detune','BH-BH Opt.','High Freq',1);
lh = legend('0) NO SRM','1a) Zero Detune, low power','1b) Zero Detune, high power', ...
    '2) NS-NS tuning','3) BH-BH tuning','location','northeast');
axis([10,3000,1e-24,1e-21]);

set(hndls(1:(end)),'LineWidth',4);
set(hndls(2),'color','r','linestyle','--')
set(gca,'fontsize',18,'fontweight','normal')
xlabel('Frequency (Hz)','fontsize',18)
ylabel('Strain noise (Hz^{-1/2})','fontsize',18)
set(lh,'fontsize',16)



return




function [y,f,score]=callGwinc(x)
shotradmode = 2;
ifo = IFOModel;
src=SourceModel;
tprmind=5;
scmode=0;
f_seis = 9;   % Seismic Wall Cutoff frequency
f_Nyq = 3000; % Fake Nyquist frequency
variate=[0,0,0,0,0,0];

   if length(x)>tprmind
     ifo.Optics.Quadrature.dc=x(6);
   end;
   switch variate(tprmind)
     case 1,
       [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x(1),x(2),x(3),x(4),x(tprmind));
     otherwise,
       [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x(1),x(2),x(3),x(4));
   end;

   y=nnn.Total;
   f=nnn.Freq;
   score = sss;
return