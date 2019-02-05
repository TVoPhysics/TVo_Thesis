% [dOpt, nLayer] = getCoatVersion(version, Ndblt, nL, nH)
%   get coating layer optical thicknesses
%     coating refractive index vector is also
%     returned if nL and nH are specified
%
% versions:
% 0 = 1/4 dH, with 1/2 wave cap
% 1 = 1/8 dH, with 1/2 wave cap
% 2 = 0.17 dH, with 1/2 wave cap

function [dOpt, nLayer] = getCoatVersion(version, Ndblt, nL, nH)

  Nlayer = 2 * Ndblt;
  dOpt = zeros(Nlayer, 1);

  % version may be a simple doublet specification
  if numel(version) == 2
    dOpt(1) = 0.5;
    dOpt(2:2:end) = version(1);
    dOpt(3:2:end) = version(2);
  elseif numel(version) == 3
    % cap, then doublet specification
    dOpt(1) = version(1);
    dOpt(2:2:end) = version(2);
    dOpt(3:2:end) = version(3);
  else
    % standard coatings, and more complex stuff
    switch version
      case 0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/4 wave
	dOpt(1) = 0.5;
	dOpt(2:end) = 0.25;

      case 1  %%%%%%%%%%%%%%%%%%%%%%%%% 1/8 - 3/8
	dOpt(1) = 0.5;
	dOpt(2:2:end) = 1 / 8;
	dOpt(3:2:end) = 3 / 8;

      case 2  %%%%%%%%%%%%%%%%%%%%%%%%% other
	dOpt(1) = 0.5;
	dOpt(2:2:end) = 0.17;
	dOpt(3:2:end) = 0.33;

      case 3  %%%%%%%%%%%%%%%%%%%%%%%%% TNI
	% 16 doublets, plus a special cap and last Ta layer (at substrate)
	%
	% doublets:
	% Si layer thickness = 250.984 nm
	% Ta layer thickness = 80.688 nm
	%
	% Si cap thickness = 29.410 nm
	% Ta 1st layer thickness = 72.677 nm
	%
	% Si refraction index = 1.46543 - i 1.8e-7
	% Ta refraction index = 2.035 - i 4e-8

	nSi = 1.446543;
	nTa = 2.035;
	lambda = 1064;
	
	dOpt(1) = 29.41 * nSi / lambda;
	dOpt(1) = 0.5;
	dOpt(2:2:end) = 80.688 * nTa / lambda;
	dOpt(3:2:end) = 250.984 * nSi / lambda;
	dOpt(end) = 72.677 * nTa / lambda;

      otherwise
	error('Unknown coating version')
    end
  end

  % set nLayer, if needed
  if nargout > 1
    nLayer = zeros(Nlayer, 1);
    nLayer(1:2:end) = nL;
    nLayer(2:2:end) = nH;
  end
  
