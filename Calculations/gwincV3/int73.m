function score = int73(f, h2, ifo, source)
% score = int73(f, h2, ifo, source)
% Takes the IFO noise and the frequency grid and make the inspiral
% score via the f^7/3 integral

 
  % try to find global f_bounce, set by gwinc in suspR
  % this is used to define the start of the integral
  global f_bounce
  if isempty(f_bounce)
    f_bounce = 9.1416;
  end

  
  MSol = ifo.Constants.MSol;
  Mpc = ifo.Constants.Mpc;
  c = ifo.Constants.c;


  g = pi*MSol/c;


  % Start the integration 0.1 Hz above the bounce mode
  n = find(f > (f_bounce+0.1));
  f = f(n);
  h2 = h2(n);


  integrand73 = (g^2) ./ ((g*f).^(7/3) .* h2);
  x73all = cumtrapz(f,integrand73);

  x73 = x73all(end);
  score.r0 = sqrt(x73*(3/20)^(5/3)*(5/(192*pi))); % Finn 96 def.
  % volume to distance conversion; from P Sutton; see T1200351
  score.r0 = (3*1.8375)^(1/3)*score.r0*MSol/Mpc;   % adjust for survey volume (in Mpc)

  % 1.4/1.4 Binary Neutron Stars -------------------------------------------------
  score.effr0ns = score.r0 * ...
    getOverlapForMassPair(source.NeutronStar, f, x73all, MSol, c);

  % Black Holes  -------------------------------------------------
  score.effr0bh = score.r0 * ...
    getOverlapForMassPair(source.BlackHole, f, x73all, MSol, c);

end

%  ------------------------------------------------- --------------------------
%
% return the score in Mpc for a inspiral source
% the source struct should contain
%   Mass1 and Mass2 in units?
%   Distance in units?
%
% x73all is the cumulative range integral

function score = getOverlapForMassPair(source, f, x73all, MSol, c)
 
  ins_mass1 = source.Mass1*MSol;
  ins_mass2 = source.Mass2*MSol;
  tot_mass = ins_mass1 + ins_mass2;
  mu = ins_mass1*ins_mass2/tot_mass;


  M0= mu^(3/5)*tot_mass^(2/5);
  distance = source.Distance;   % in Mpc
  z = distance*75000/c;
  MM = M0*(1+z);


  f2max = 2 * (1980/(z+1)) * (MSol/tot_mass);
  if (f2max < f(1))
    error('f_isco is less than f_low')
  end

  x73 = x73all(end);
  zetafmax = interp1(f, x73all, f2max) / x73;

  % debug: make sure point lies on curve
  %semilogx(f, x73all, f2max, x73 * zetafmax, 'x');
  %pause
  
  score = sqrt(zetafmax) * (MM/1.2188/MSol)^(5/6);

end
