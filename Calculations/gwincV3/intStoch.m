function score = intStoch(f, h2, OverlapMode, ifo, source)
%
% Takes the IFO noise and the frequency grid and make the stochastic
% score via the f^7/3 integral

npow = source.Stochastic.powerlaw;

% try to find global f_bounce, set by gwinc in suspR
% this is used to define the start of the integral
global f_bounce
if isempty(f_bounce)
  f_bounce = 9.1416;
end

% Start the integration 0.1 Hz above the bounce mode
n = find(f > (f_bounce+0.1));
f = f(n);
h2 = h2(n);

% use global nse for overlap if OverlapMode is not 0
global nse;
switch OverlapMode
  case 0
   nse.ovlp = ovlp(f, ifo);
   integrandSt = nse.ovlp.^2.*f.^(npow-6)./(h2.^2);
  otherwise
    integrandSt = nse.ovlp.^2.*f.^(npow-6)./(h2.^2);
end;
x = trapz(f,integrandSt);

x = 2*x;				                      % positive & negative frequencies
H0 = source.Stochastic.Hubble/(ifo.Constants.Mpc/10^3);       % H0 in units Hz
time = source.Stochastic.integration_time * ifo.Constants.yr; % integration time
x = x*(3*H0^2/(10*pi^2))^2*time;	                      % 

% Omega*h_100^2
% Factor 2.5 is for 99% unified confidence interval on Gaussian test statistic. 
score = 2.5/sqrt(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = ovlp(f, ifo)
% OVLP - "overlap reduction function" for stochastic signal
% searches in LIGO. Reference Allen Les Houche lectures.   

ckm = ifo.Constants.c/10^3;		% km/s
d = 3001;                 		% separation between LIGO sites (km)
a = 2*pi*f*d/ckm;

j0 = sin(a)./a;
j1 = j0./a - cos(a)./a;
j2 = 3*j1./a - j0;

j1 = j1./a;
j2 = j2./a.^2;

% Coefficients from Allen Les Houche lectures correspond to LIGO
% detector relative arm orientations
g = -0.124842*j0 - 2.90014*j1 + 3.00837*j2;  

return
