function n = subbrownian(f, ifo)
% SUBBROWNIAN - strain noise psd arising from the Brownian thermal noise 
% due to mechanical loss in the substrate material

  wITM = ifo.Optics.ITM.BeamRadius;
  wETM = ifo.Optics.ETM.BeamRadius;
  Y = ifo.Materials.Substrate.MirrorY;
  sigma = ifo.Materials.Substrate.MirrorSigma;

  c2 = ifo.Materials.Substrate.c2;
  n = ifo.Materials.Substrate.MechanicalLossExponent;
  alphas = ifo.Materials.Substrate.Alphas;
  L = ifo.Infrastructure.Length;
  kBT = ifo.Constants.kB*ifo.Constants.Temp;

  % Bulk substrate contribution
  phibulk = c2.*f.^n;

  [cITM, aITM] = subbrownianFiniteCorr(ifo, 'ITM');
  [cETM, aETM] = subbrownianFiniteCorr(ifo, 'ETM');
  cbulk = 8 * kBT * (aITM + aETM) .* phibulk ./ (2 * pi * f);

  % Surface loss contribution
  % csurfETM = alphas/(Y*pi*wETM^2);
  % csurfITM = alphas/(Y*pi*wITM^2);

  csurfETM = alphas*(1-2*sigma)/((1-sigma)*Y*pi*wETM^2);
  csurfITM = alphas*(1-2*sigma)/((1-sigma)*Y*pi*wITM^2);
  csurf = 8 * kBT * (csurfITM + csurfETM) ./ (2 * pi * f);

  % account for 2 ITM and 2 ETM, and convert to strain whith 1/L^2
  n = 2 * (csurf + cbulk) / L^2;

