function n = subtherm(f, ifo)
% SUBTHERM - noise psd arising from thermoelastic fluctuations in mirror

  wITM = ifo.Optics.ITM.BeamRadius;
  wETM = ifo.Optics.ETM.BeamRadius;
  sigma = ifo.Materials.Substrate.MirrorSigma;

  L = ifo.Infrastructure.Length;
  kBT = ifo.Constants.kB*ifo.Constants.Temp;

  rho = ifo.Materials.Substrate.MassDensity;	%
  kappa = ifo.Materials.Substrate.MassKappa;	% thermal conductivity
  alpha = ifo.Materials.Substrate.MassAlpha;	% thermal expansion
  CM = ifo.Materials.Substrate.MassCM;		% heat capacity @ constant mass
  Temp = ifo.Constants.Temp;		% temperature

  S0 = 8*(1+sigma)^2*kappa*alpha^2*Temp*kBT;	% note kBT has factor Temp
  S0 = S0/(sqrt(2*pi)*(CM*rho)^2);
  SITM = S0/(wITM/sqrt(2))^3;		% LT 18 less factor 1/omega^2
  SETM = S0/(wETM/sqrt(2))^3;		% LT 18 less factor 1/omega^2

  % Corrections for finite test masses:
  SITM = SITM * subthermFiniteCorr(ifo, 'ITM');
  SETM = SETM * subthermFiniteCorr(ifo, 'ETM');
  
  % 2 mirrors each type, factor omega^2, dimensionless strain
  n = 2*(SITM + SETM)./(2*pi*f*L).^2;
