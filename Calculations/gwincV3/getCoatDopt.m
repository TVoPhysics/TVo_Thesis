% dOpt = getCoatDopt(ifo, opticName)
% dOpt = getCoatDopt(ifo, T, dL, dCap)
% dOpt = getCoatDopt([nS, nL, nH], T, dL, dCap)
%   get coating layer optical thicknesses to match desired transmission
%
% ifo = gwinc IFO model, or vector of 3 refractive indexes [nS, nL, nH]
%   nS, nL, nH = refractive index of substrate, low and high-n materials
%     nL is used for the even layers, nH for odd layers
%     this algorithm may fail if nS > nL
% opticName = name of the Optic struct to use for T, dL and dCap
% T = power transmission of coating
% dL = optical thickness of low-n layers
%   high-n layers have dH = 0.5 - dL
% dCap = first layer (low-n) thickness (default 0.5)
%
% dOpt = optical thickness vector Nlayer x 1
%
% Examples:
% dOpt = getCoatDopt(ifo, 'ITM');
% dOpt = getCoatDopt(ifo, 0.014, 0.3);
% dOpt = getCoatDopt(ifo, 0.014, 0.3, 0.5);
% dOpt = getCoatDopt([1.45, 1.45, 2.06], 0.014, 0.3, 0.5);

function dOpt = getCoatDopt(ifoArg, T, dL, dCap)

  % check arguments
  if nargin < 3
    % T has been used at the name of an optic
    opt = ifoArg.Optics.(T);
    
    T = opt.Transmittance;
    dL = opt.CoatingThicknessLown;
    dCap = opt.CoatingThicknessCap;
  elseif nargin < 4
    dCap = 0.5;
  end
  
  % get IFO model stuff (or create it for other functions)
  if numel(ifoArg) == 3
    nS = ifoArg(1);
    nL = ifoArg(2);
    nH = ifoArg(3);

    ifo.Materials.Substrate.RefractiveIndex = nS;
    ifo.Materials.Coating.Indexlown = nL;
    ifo.Materials.Coating.Indexhighn = nH;
  else
    ifo = ifoArg;
    pS = ifo.Materials.Substrate;
    pC = ifo.Materials.Coating;
  
    nS = pS.RefractiveIndex;
    nL = pC.Indexlown;  
    nH = pC.Indexhighn;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % find number of quarter-wave layers required, as first guess
  nR = nH / nL;
  a1 = (2 - T + 2 * sqrt(1 - T)) / (nR * nH * T);
  Ndblt = ceil(log(a1) / (2 * log(nR)));
  
  % search through number of doublets to find Ndblt
  % which gives T lower than required
  dH = 0.5 - dL;
  Tn = getTrans(ifo, Ndblt, dL, dH, dCap, dH);
  while Tn < T && Ndblt > 1
    % strange, but T is too low... remove doublets
    Ndblt = Ndblt - 1;
    Tn = getTrans(ifo, Ndblt, dL, dH, dCap, dH);
  end
  while Tn > T && Ndblt < 1e3
    % add doublets until T > tN
    Ndblt = Ndblt + 1;
    Tn = getTrans(ifo, Ndblt, dL, dH, dCap, dH);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % tweak bottom layer
  delta = 0.01;
  dScan = (0:delta:0.25)';
  dTweak = getTweak(ifo, T, Ndblt, dL, dH, dCap, dScan, 5);

  if isempty(dTweak)
    if nS > nL
      error('Coating tweak layer not sufficient since nS > nL.')
    else
      error('Coating tweak layer not found... very strange.')
    end
  end
  
  % now that we are close, get a better result with a linear fit
  delta = 0.001;
  dScan = (dTweak - 3 * delta:delta:dTweak + 3 * delta)';
  [dTweak, Td] = getTweak(ifo, T, Ndblt, dL, dH, dCap, dScan, 3);
  
  % negative values are bad
  if dTweak < 0.01
    dTweak = 0.01;
  end
  
  % check the result
  if abs(log(Td / T)) > 1e-3
    warning('Exact coating tweak layer not found... %g%% error.', ...
      abs(log(Td / T)))
  end

  %%%%%%%%%%%%%%%%%%%%%%%%
  % return dOpt vector
  dOpt = zeros(2 * Ndblt, 1);
  dOpt(1) = dCap;
  dOpt(2:2:end) = dH;
  dOpt(3:2:end) = dL;
  dOpt(end) = dTweak;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = getTrans(ifo, Ndblt, dL, dH, dCap, dTweak)
  
  % the optical thickness vector
  dOpt = zeros(2 * Ndblt, 1);
  dOpt(1) = dCap;
  dOpt(2:2:end) = dH;
  dOpt(3:2:end) = dL;
  
  N = numel(dTweak);
  T = zeros(N, 1);
  for n = 1:N
    dOpt(end) = dTweak(n);
    r = getCoatRefl(ifo, dOpt);
    T(n) = 1 - abs(r^2);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dTweak, Td] = getTweak(ifo, T, Ndblt, dL, dH, dCap, dScan, Nfit)

  % tweak bottom layer
  Tn = getTrans(ifo, Ndblt, dL, dH, dCap, dScan);
  pf = polyfit(dScan, Tn - T, Nfit);
  rts = roots(pf);
  dTweak = min(rts(imag(rts) == 0 & rts > 0));

  if isempty(dTweak)
    Td = 0;
    return
  end
  
  % compute T for this dTweak
  Td = getTrans(ifo, Ndblt, dL, dH, dCap, dTweak);

  % plot for debugging
%   plot(dScan, [Tn - T, polyval(pf, dScan)], dTweak, Td - T, 'ro')
%   grid on
%   legend('T exact', 'T fit', 'T at dTweak')
%   title(sprintf('%d doublets', Ndblt))
%   pause(1)
end
