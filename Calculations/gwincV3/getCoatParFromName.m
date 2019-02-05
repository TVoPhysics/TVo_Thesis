% [wBeam, dOpt] = getCoatParFromName(ifo, optName)
%   return parameters relevant to coating calculations given an optic's name

function [wBeam, dOpt] = getCoatParFromName(ifo, optName)

  if ~isfield(ifo.Optics, optName)
    error('Unknown optic specified.  optName = %s', optName);
  end
  opt = ifo.Optics.(optName);

  % beam size
  wBeam = opt.BeamRadius;

  % coating layers
  if isfield(opt, 'CoatLayerOpticalThickness')
    dOpt = opt.CoatLayerOpticalThickness;
  else
    dOpt = getCoatDopt(ifo, optName);
  end
  
