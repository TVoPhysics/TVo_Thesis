%----------------------------------------------------------------
% function [x,y] = VarySR3_95MM_wSqz_MMTTon_maxtem4(noplot)
% Matlab function to plot Finesse output data
% Usage: 
%   [x,y] = VarySR3_95MM_wSqz_MMTTon_maxtem4    : plots and returns the data
%   [x,y] = VarySR3_95MM_wSqz_MMTTon_maxtem4(1) : just returns the data
%           VarySR3_95MM_wSqz_MMTTon_maxtem4    : just plots the data
% Created automatically Tue Sep 13 22:09:04 2016
% by Finesse 2.1 (2.1-2-gb87c29b), 16.05.2016
%----------------------------------------------------------------
function [x,y] = VarySR3_95MM_wSqz_MMTTon_maxtem4(noplot)

data = load('VarySR3_95MM_wSqz_MMTTon_maxtem4.out');
[rows,cols]=size(data);
x=data(:,1);
y=data(:,2:cols);
mytitle='VarySR3\_95MM\_wSqz\_MMTTon\_maxtem4                Tue Sep 13 22:09:04 2016';
if (nargin==0)

figure('name','VarySR3_95MM_wSqz_MMTTon_maxtem4');

h1=subplot(2,1,1);
plot(x, y(:,1), x, y(:,3), x, y(:,5), x, y(:,7), x, y(:,9), x, y(:,11), x, y(:,13), x, y(:,15), x, y(:,17), x, y(:,19), x, y(:,21), x, y(:,23), x, y(:,25), x, y(:,27), x, y(:,29), x, y(:,31));
legend('NSR\_with\_RP nOMC\_AROC\_trans', 'nSRMHRaTEM00 nSRMHRa', 'nSRMHRaTEM01 nSRMHRa', 'nSRMHRaTEM02 nSRMHRa', 'SRCoutx nIBAin', 'SRCouty nIBAin', 'SRMYqx nSRMHRa', 'SRMYqy nSRMHRa', 'ITMXqx nITMX2', 'ITMXqy nITMX2', 'ITMYqx nITMY2', 'ITMYqy nITMY2', 'OMCqx nOMC\_HROC\_refl', 'OMCqy nOMC\_HROC\_refl', 'OFIqx nIBAin', 'OFIqy nIBAin');
set(gca, 'YScale', 'lin');
ylabel('Re ');
set(gca, 'XLim', [5 5000]);
set(gca, 'XScale', 'log');
xlabel('f [Hz] (darm)');
grid on;

h2=subplot(2,1,2);
plot(x, y(:,2), x, y(:,4), x, y(:,6), x, y(:,8), x, y(:,10), x, y(:,12), x, y(:,14), x, y(:,16), x, y(:,18), x, y(:,20), x, y(:,22), x, y(:,24), x, y(:,26), x, y(:,28), x, y(:,30), x, y(:,32));
legend('NSR\_with\_RP nOMC\_AROC\_trans', 'nSRMHRaTEM00 nSRMHRa', 'nSRMHRaTEM01 nSRMHRa', 'nSRMHRaTEM02 nSRMHRa', 'SRCoutx nIBAin', 'SRCouty nIBAin', 'SRMYqx nSRMHRa', 'SRMYqy nSRMHRa', 'ITMXqx nITMX2', 'ITMXqy nITMX2', 'ITMYqx nITMY2', 'ITMYqy nITMY2', 'OMCqx nOMC\_HROC\_refl', 'OMCqy nOMC\_HROC\_refl', 'OFIqx nIBAin', 'OFIqy nIBAin');
ylabel('Im ');
set(gca, 'YScale', 'lin');
set(gca, 'XLim', [5 5000]);
set(gca, 'XScale', 'log');
xlabel('f [Hz] (darm)');
grid on;

subplot(2,1,1);
title(mytitle);
end

switch nargout
 case {0}
  clear x y;
 case {2}
 otherwise
  error('wrong number of outputs');
end
