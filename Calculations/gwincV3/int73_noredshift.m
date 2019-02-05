function score = int73(f, h2, ifo, Mass)
%
% Takes the IFO noise and the frequency grid and make the inspiral
% score via the f^7/3 integral, without correcting for redshift

MSol = ifo.Constants.MSol;
Mpc = ifo.Constants.Mpc;
c = ifo.Constants.c;
%fInspiralMin = ifo.Constants.fInspiralMin;
fInspiralMin = 9.1416; % Hz   arbitrary choice I guess

g = pi*MSol/c;

% Start the integration 0.1 Hz above the bounce mode
n = find(f > fInspiralMin);
f = f(n);
h2 = h2(n);

integrand73 = (g^2) ./ ((g*f).^(7/3) .* h2);

%%  1.4/1.4 Binary Neutron Stars -------------------------------------------------
ins_mass1 = Mass*MSol;
ins_mass2 = Mass*MSol;
tot_mass = ins_mass1 + ins_mass2;
mu = ins_mass1*ins_mass2/tot_mass;

M0= mu^(3/5)*tot_mass^(2/5);
distance = 0;   % in Mpc
z = distance*75000/c;
MM = M0*(1+z);

f2max = 2 * (1980/(z+1)) * (MSol/tot_mass);

x73all = cumtrapz(f,integrand73);
x73 = x73all(end);

nn = find(f < f2max);
zetafmax_num = x73all(nn(end));

if (f2max < f(1))
  error('f_isco is less than f_low')
end

zetafmax = zetafmax_num/x73;

score.r0 = sqrt(x73*(3/20)^(5/3)*(5/(192*pi))); % Finn 96 def.
score.r0 = (3*1.8375)^(1/3)*score.r0*MSol/Mpc;     % adjust for survey volume (in Mpc)  ...
%correct averageing according to Patrick Sutton's note (Daniel)
score.effr0 = 2.2643.* score.r0 * sqrt(zetafmax) * (MM/1.2188/MSol)^(5/6); %2.2643 is to convert to horizont dist
%  ------------------------------------------------- --------------------------



%  ------------------------------------------------- --------------------------
