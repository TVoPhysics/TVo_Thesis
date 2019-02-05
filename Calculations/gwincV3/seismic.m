function n = seismic(f,ifo)
% SEISMIC - seismic noise psd at frequencies f for given ifo.
%
% Modified to include realistic SEI + SUS models (Rana, 11/2005)

  % Interpolate the log10 onto the ifo frequency array
  n = interp1(ifo.Seismic.darmseis_f, ...
    log10(ifo.Seismic.darmseis_x), f, 'cubic', -30);

  % Convert into Strain
  n = 10.^(n) / ifo.Infrastructure.Length;

  % Convert to PSD from ASD
  n = n.^2;
