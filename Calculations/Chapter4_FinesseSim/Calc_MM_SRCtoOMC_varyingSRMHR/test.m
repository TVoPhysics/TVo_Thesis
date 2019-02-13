%----------------------------------------------------------------
% function [x,y] = test(noplot)
% Matlab function to plot Finesse output data
% Usage: 
%   [x,y] = test    : plots and returns the data
%   [x,y] = test(1) : just returns the data
%           test    : just plots the data
% Created automatically Wed Sep  7 13:42:53 2016
% by Finesse 2.1 (2.1-2-gb87c29b), 16.05.2016
%----------------------------------------------------------------
function [x,y] = test(noplot)

data = load('test.out');
[rows,cols]=size(data);
x=data(:,1);
y=data(:,2:cols);
mytitle='test                Wed Sep  7 13:42:53 2016';
if (nargin==0)

figure('name','test');

h1=subplot(2,1,1);
plot(x, y(:,1));
legend('NSR\_with\_RP nOMC\_AROC\_trans');
set(gca, 'YScale', 'lin');
ylabel('Re ');
set(gca, 'XLim', [5 5000]);
set(gca, 'XScale', 'log');
xlabel('f [Hz] (darm)');
grid on;

h2=subplot(2,1,2);
plot(x, y(:,2));
legend('NSR\_with\_RP nOMC\_AROC\_trans');
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
