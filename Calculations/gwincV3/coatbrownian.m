function n = coatbrownian(f,ifo)
% COAT - effect of optical coating on Brownian thermal noise
%   returns strain noise power spectrum in 1 / Hz
%
% Added by G Harry 8/3/02 from work by Nakagawa, Gretarsson, et al.
% Expanded to reduce approximation, GMH 8/03
% Modified to return strain noise, PF 4/07
% Modified to accept coating with non-quarter-wave layers, mevans 25 Apr 2008

  % Constants
  L = ifo.Infrastructure.Length;

  % compute noise power from one ITM and one ETM
  SbrITM  = getCoatBrownian(f, ifo, 'ITM');
  SbrETM  = getCoatBrownian(f, ifo, 'ETM');

  n = 2 * (SbrITM + SbrETM) / L^2;
