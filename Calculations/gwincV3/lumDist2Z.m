function [z,dC,dM,dA,dL]=lumDist2Z(dL,H0,OmegaM,OmegaLambda,OmegaK)
% returns
%  redshift                     z
%  Comoving distance            dC
%  Transverse comoving distance dM
%  Angular diameter distance    dA
%  Luminosity distance          dL
%
% all distances are in Gpc
% Luminosity Distance dL is required, all other arguments are optional
%
% Stefan Ballmer 2014/8/1
% sballmer@syr.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference data (overwritten by arguments)
ref.H0    =67.9;   % ref:  arXiv:1303.5076
ref.OmegaM=0.307;   % ref:  arXiv:1303.5076

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z values for numerical integration
dz0=1e-3;
z0=0:dz0:10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% argument checking
if ~exist('dL','var')
    error('Luminosity Distance dL is required argument');
end
if ~exist('H0','var')
    H0=ref.H0;
end
if ~exist('OmegaM','var')
    OmegaM=ref.OmegaM;
end
if ~exist('OmegaLambda','var')
    if ~exist('OmegaK','var')
%        disp('Flat universe assumed');
        OmegaLambda=1-OmegaM;
        OmegaK=0;
    else
        OmegaLambda=1-OmegaM-OmegaK;
    end
else
    if ~exist('OmegaK','var')
        OmegaK=1-OmegaM-OmegaLambda;
    elseif abs(OmegaM+OmegaK+OmegaLambda-1)>1e-6
        error('Total of energies must be =1');
    end
end
if OmegaK~=0
    disp('Universe is not flat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed constants
Mpc=3.08567758e22; % meter
Gpc=1e3*Mpc; % meter
km=1e3;            % meter
H0SI=H0*km/Mpc;    % Hubble constant in SI units
c=299792458;       % speed of light
dHSI=c/H0SI;       % Hubble distance in meters
dH=dHSI/Gpc;     % Hubble distance in Gpc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical integration
E0=sqrt(OmegaM.*(1+z0).^3+OmegaK.*(1+z0).^2+OmegaLambda);
I0=cumtrapz(1./E0)*dz0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the results

%Comoving distance
dC0=dH*I0;

%Transverse comoving distance
if OmegaK==0
    dM0=dC0;
elseif OmegaK>0
    g=sqrt(OmegaK);
    dM0=(dH/g).*sinh(g.*I0);
else
    g=sqrt(abs(OmegaK));
    dM0=(dH/g).*sin(g.*I0);
end

%Angular diameter distance
dA0=dM0./(1+z0);

%Luminosity distance
dL0=dM0.*(1+z0);

% red shift and other results
z=interp1(dL0,z0,dL);
dC=interp1(dL0,dC0,dL);
dM=interp1(dL0,dM0,dL);
dA=interp1(dL0,dA0,dL);






