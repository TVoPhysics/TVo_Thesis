function n = thermooptic(f, ifo)
% COATTHERMO - effect of thermoelastic and thermorefractive fluctuations in
% the mirror dielectric coating
%   returns strain noise power spectrum in 1 / Hz
%
% Added by G Harry 8/27/03 from work by Fejer, Rowan, Braginsky, et al
% thermoelastic and thermorefractive effects combined coherently in 2006
%
% Reduced approximations and combined TE and TR
%   effects with correct sign, mevans 25 Apr 2008


  % Constants
  L = ifo.Infrastructure.Length;

  % compute noise power from one ITM and one ETM
  StoITM  = getCoatThermoOptic(f, ifo, 'ITM');
  StoETM  = getCoatThermoOptic(f, ifo, 'ETM');

  n = 2 * (StoITM + StoETM) / L^2;
