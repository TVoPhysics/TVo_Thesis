% ifo = precompIFO(ifo, PRfixed)
%   add precomputed data to the IFO model
%
% To prevent recomputation of these precomputed data, if the
% ifo argument contains ifo.gwinc.PRfixed, and this matches
% the argument PRfixed, no changes are made.
%
% (mevans June 2008)

function ifo = precompIFO(ifo, PRfixed)
  
  % check PRfixed
  if isfield(ifo, 'gwinc')
    % && isfield(ifo.gwinc, 'PRfixed') && ifo.gwinc.PRfixed == PRfixed
    return
  end
  ifo.gwinc.PRfixed = PRfixed;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DERIVED OPTICS VALES
  % Calculate optics' parameters
  ifo.Materials.MirrorMass = ...
    pi*ifo.Materials.MassRadius^2*ifo.Materials.MassThickness;
  ifo.Materials.MirrorMass = ifo.Materials.MirrorMass* ...
    ifo.Materials.Substrate.MassDensity;		% Kg
  ifo.Optics.ITM.Thickness = ifo.Materials.MassThickness;

  % coating layer optical thicknesses - mevans 2 May 2008
  ifo.Optics.ITM.CoatLayerOpticalThickness = getCoatDopt(ifo, 'ITM');
  ifo.Optics.ETM.CoatLayerOpticalThickness = getCoatDopt(ifo, 'ETM');

  % compute power on BS
  [pbs, finesse, prfactor, Tpr] = BSPower(ifo, PRfixed);
  ifo.gwinc.pbs = pbs;
  ifo.gwinc.finesse = finesse;
  ifo.gwinc.prfactor = prfactor;
  ifo.Optics.PRM.Transmittance = Tpr;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SAVED DATA  
  % precompute bessels zeros for finite mirror corrections
  global besselzeros;
  if isempty(besselzeros)
    % load saved values, or just compute them
    try
      load besselzeros
    catch
      besselzeros = besselzero(1, 300, 1);
    end
  end
  ifo.Constants.BesselZeros = besselzeros;

  % Seismic noise term is saved in a .mat file defined in your respective IFOModel.m
  % It is loaded here and put into the ifo structure.
  load(ifo.Seismic.darmSeiSusFile)

  ifo.Seismic.darmseis_f = darmseis_f;
  ifo.Seismic.darmseis_x = darmseis_x;

