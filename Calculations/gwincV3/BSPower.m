function [pbs, finesse, prfactor, Tpr] = BSPower(ifo, PRfixed)
% BSPower - Find and return power on ifo beamsplitter limited by 
% thermal lensing effects, finesse and power recycling factor of ifo
% [pbs,finesse,prfactor] = BSPower(ifo)

  % Modified to use a fixed PR factor   Rana, Feb 07

  % constants
  c       = ifo.Constants.c;
  pin     = ifo.Laser.Power;
  lambda  = ifo.Laser.Wavelength;
  t1      = sqrt(ifo.Optics.ITM.Transmittance);
  r1      = sqrt(1 - ifo.Optics.ITM.Transmittance);
  t2      = sqrt(ifo.Optics.ETM.Transmittance);
  r2      = sqrt(1 - ifo.Optics.ETM.Transmittance);
  t3      = sqrt(ifo.Optics.SRM.Transmittance);
  t5      = sqrt(ifo.Optics.PRM.Transmittance);
  r5      = sqrt(1 - ifo.Optics.PRM.Transmittance);
  wl      = 2*pi*c/lambda;
  lrec    = ifo.Optics.SRM.CavityLength;
  effic   = ifo.Optics.PhotoDetectorEfficiency;
  loss    = ifo.Optics.Loss;                          % single TM loss
  bsloss  = ifo.Optics.BSLoss;
  acoat   = ifo.Optics.ITM.CoatingAbsorption;
  pcrit   = ifo.Optics.pcrit;

  % Finesse, effective number of bounces in cavity, power recycling factor
  finesse = 2*pi/(t1.^2 + 2*loss);        % arm cavity finesse
  neff    = 2*finesse/pi;

  % Arm cavity reflectivity with finite loss
  rarm = r1 - (t1^2) * r2 * sqrt(1-2*loss) / (1 - r1*r2*sqrt(1-2*loss));

  if (PRfixed == 1)
    Tpr = ifo.Optics.PRM.Transmittance;  % use given value
  else
    %prfactor = 1/(2*loss * neff + bsloss);         % power recycling factor
    Tpr = 1-(rarm*sqrt(1-bsloss))^2; % optimal recycling mirror transmission
    t5 = sqrt(Tpr);
    r5 = sqrt(1 - Tpr);
  end
  prfactor = t5^2 / (1 + r5 * rarm * sqrt(1-bsloss))^2;

  pbs = pin*prfactor;			% bs power from input power

  asub = 1.3*2*ifo.Optics.ITM.Thickness*ifo.Optics.SubstrateAbsorption;
  pbsl = 2*pcrit/(asub+acoat*neff);	% bs power limited by thermal lensing

  %pbs = min([pbs,pbsl]);
  if pbs > pbsl
    disp('P_BS exceeds BS Thermal limit!')
  end
