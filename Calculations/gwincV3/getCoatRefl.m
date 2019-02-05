% rCoat = getCoatRefl(ifo, dOpt)
%   returns amplitude reflectivity, with phase, of a coating
%
% ifo = parameter struct from IFOmodel.m
% dOpt   = optical thickness / lambda of each layer
%        = geometrical thickness * refractive index / lambda
%
% This code is taken from ewa/multidiel1.m

function rCoat = getCoatRefl(ifo, dOpt)

  pS = ifo.Materials.Substrate;
  pC = ifo.Materials.Coating;
  
  nS = pS.RefractiveIndex;
  nL = pC.Indexlown;  
  nH = pC.Indexhighn;
  
  Nlayer = numel(dOpt);

  % refractive index of input, coating, and output materials
  nAll = zeros(Nlayer + 2, 1);
  nAll(1) = 1;  % vacuum input
  nAll(2:2:end) = nL;
  nAll(3:2:end) = nH;
  nAll(end) = nS; % substrate output

  % reflectivity of each interface
  rInterface = (nAll(1:end-1) - nAll(2:end)) ./ (nAll(1:end-1) + nAll(2:end));

  % combine layers as small FP cavities, starting with last reflectivity
  rCoat = rInterface(end);
  for n = numel(dOpt):-1:1
    % one-way phase in this layer
    phi = 2 * pi * dOpt(n);

    % accumulate reflecitivity
    rCoat = rCoat .* exp(-2 * i * phi);
    rCoat = (rInterface(n) + rCoat) ./ (1 + rInterface(n) * rCoat);
  end
