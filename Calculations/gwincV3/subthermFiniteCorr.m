function coeff = subthermFiniteCorr(ifo, opticName)
% CFTM20 - finite size test mass correction to noise amplitude coefficient
% (Liu & Thorne gr-qc/0002055 equation 46)
% 
% Equation references to Bondu, et al. Physics Letters A 246 (1998)
% 227-236 (hereafter BHV) or Liu and Thorne gr-qc/0002055 (hereafter LT)

  % extract some numbers
  a = ifo.Materials.MassRadius;
  h = ifo.Materials.MassThickness;
  w = ifo.Optics.(opticName).BeamRadius;
  sigma = ifo.Materials.Substrate.MirrorSigma;
  zeta = ifo.Constants.BesselZeros;

  % do the work
  j0m = besselj(0,zeta);
  r0 = w/sqrt(2);			% LT uses power e-folding
  km = zeta/a;

  Qm = exp(-2*km*h);			% LT eq. 35a

  pm = exp(-(km*r0).^2/4)./(pi*(a*j0m).^2);	% LT 37

  c0 = 6*(a/h)^2*sum(j0m.*pm./zeta.^2);	% LT 32
  c1 = -2*c0/h;				% LT 32
  p0 = 1/(pi*a^2);			% LT 28
  c1 = c1 + p0/(2*h);			% LT 40

  coeff = (1-Qm).*((1-Qm).*(1+Qm)+8*h*km.*Qm);
  coeff = coeff + 4*(h*km).^2.*Qm.*(1+Qm);
  coeff = coeff.*km.*(pm.*j0m).^2.*(1-Qm);
  coeff = coeff./((1-Qm).^2-4*(h*km).^2.*Qm).^2;
  coeff = sum(coeff) + h*c1^2/(1+sigma)^2;
  coeff = coeff*(sqrt(2*pi)*r0)^3*a^2;	% LT 46
