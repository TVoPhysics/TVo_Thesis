function varargout = gwinc(flo,fhi,ifoin,sourcein,varargin)
% GWINC   Calculates strain noise due to various noise sources, for a
% specified set of interferometer parameters. Also evaluates the
% sensitivity of the interferometer to the detection of several potential 
% gravitational wave sources. Usage:
%
%      VARARGOUT = GWINC(FLO,FHI,IFO,SOURCE,VARARGIN)
%
%      FLO, FHI = minimum and maximum frequencies between which
%                  calculations are made
%      IFO       = structure containing interferometer parameters
%      SOURCE    = structure containing source parameters
%
% Optional input arguments (the last 4 override IFO parameters):
%      VARARGIN{1}: PLOT_FLAG set to 4 for score, only calculating shotrad
%                                    3 for score and plots
%                                    2 for score only
%                                    1 to make plots but no score
%                                    else 0 (DEF)
%      VARARGIN{2}: LASER POWER -> ifo.Laser.Power
%      VARARGIN{3}: SRC PHASE   -> ifo.Optics.SRM.Tunephase
%      VARARGIN{4}: SRM TRANS   -> ifo.Optics.SRM.Transmittance
%      VARARGIN{5}: ITM TRANS   -> ifo.Optics.ITM.Transmittance
%      VARARGIN{6}: PRM TRANS   -> ifo.Optics.PRM.Transmittance
%      VARARGIN{7}: HOMO PHASE  -> ifo.Optics.Quadrature.dc
%
% Optional output arguments
%      VARARGOUT{1}: SCORE  structure containing source sensitivities
%      VARARGOUT{2}: NOISE  structure containing noise terms
%
% Ex.1    [score,noise] = gwinc(5,5000,IFOModel,SourceModel,1)



if (nargin < 4)
  error('usage: gwinc(flo,fhi,ifo,source,...);');
end

% avoid modifying arguments
ifo = ifoin;
source = sourcein;

% -------------------------------------------------------
% parse arguments

fig = 0;
makescore = 0;
modeSR = 0;
% PRfixed: set to 1 for a fixed PRM transmission; 0 to allow optimization,
% which occurs in BSPower.m
PRfixed = 1;   

% Parse varargin to decide to make plots or scores or both or neither
if (nargin > 4)
  if varargin{1} > 1
    makescore = 1;
  end
  if (varargin{1} == 1 | varargin{1} == 3)
    fig = 1;
  end
  if varargin{1} == 4
    modeSR = 1;
  end
end

% Adjust these parameters as command line arguments
if nargin > 5
   ifo.Laser.Power = varargin{2};
end
if nargin > 6
  ifo.Optics.SRM.Tunephase = varargin{3};
end
if nargin > 7
  ifo.Optics.SRM.Transmittance  = varargin{4};
end
if nargin > 8
  ifo.Optics.ITM.Transmittance  = varargin{5};
end
if nargin > 9
  PRfixed = 1;
  ifo.Optics.PRM.Transmittance  = varargin{6};
end
if nargin > 10
    ifo.Optics.Quadrature.dc = varargin{7};
end
if nargin > 11
  error('Too many arguments to gwinc')
end

% --------------------------------------------------------
% add some precomputed info to the ifo struct
ifo = precompIFO(ifo, PRfixed);

pbs = ifo.gwinc.pbs;
finesse = ifo.gwinc.finesse;
prfactor = ifo.gwinc.prfactor;
armpower = finesse*(2/pi)*(pbs/2)*(1/1000); % kW
PowAbsITM = [finesse*(2/pi)*ifo.Optics.ITM.CoatingAbsorption/2; ...
             ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption/2]*pbs;
ifo.TCS.M = [ifo.TCS.s_cc ifo.TCS.s_cs;...
             ifo.TCS.s_cs ifo.TCS.s_ss];         
S_uncorr = transpose(PowAbsITM)*ifo.TCS.M*PowAbsITM;         
TCSeff = 1 - sqrt(ifo.TCS.SRCloss/S_uncorr);
thermalLoad.ITM = sum(PowAbsITM);
thermalLoad.BS = ifo.Materials.MassThickness * ifo.Optics.SubstrateAbsorption * pbs;

if (ifo.Laser.Power*prfactor ~= pbs)
  disp(sprintf('Warning: lensing limits input power to %7.2f W',...
  		pbs/prfactor));
end

% --------------------------------------------------------
% Compute all noise sources and sum for plotting and pulsar benchmark

% Frequency grid on which everything is calculated
f = logspace(log10(flo),log10(fhi),3000);

% Physical constants
c = ifo.Constants.c;
MSol = ifo.Constants.MSol;
Mpc = ifo.Constants.Mpc;
yr = ifo.Constants.yr;

% global noise source struct can be used with plot option = 4
%   modeSR = 0 => compute all noises and store them in nse
%   modeSR other => compute only quantum noise, use old nse for others
global nse

% Compute noises
y1 = shotrad(f,ifo);
switch modeSR
case 0
 y3 = suspR(f,ifo);
 y4 = gas(f,ifo);
 y5 = subbrownian(f,ifo);                      % substrate Brownian
 y6 = coatbrownian(f,ifo);                     % Coating Brownian
 y8 = subtherm(f,ifo);                         % substrate thermo-elastic
 y9 = gravg(f,ifo);
 y10 = seismic(f,ifo);
 y11 = thermooptic(f,ifo);                     % coating thermo-optic (TE + TR)
 y2 = y5 + y6 + y8 + y11;                      % total mirror thermal 
otherwise
 y3 = nse.SuspThermal;
 y4 = nse.ResGas;
 y5 = nse.MirrorThermal.SubBrown;              % substrate Brownian
 y6 = nse.MirrorThermal.CoatBrown;             % Coating Brownian
 y8 = nse.MirrorThermal.SubTE;                 % substrate thermo-elastic
 y9 = nse.Newtonian;
 y10 = nse.Seismic;
 y11 = nse.MirrorThermal.CoatTO;               % coating thermo-optic
 y2 = nse.MirrorThermal.Total;                 % total mirror thermal 
end
ys = y1 + y2 + y3 + y4 + y9 + y10;             % sum of noise sources

% --------------------------------------------------------
% output graphics

if (fig ~= 0)
  % Report input parameters
  fprintf('Laser Power:            %7.2f Watt\n',ifo.Laser.Power);
  fprintf('SRM Detuning:           %7.2f degree\n',ifo.Optics.SRM.Tunephase*180/pi);
  fprintf('SRM transmission:       %9.4f\n',ifo.Optics.SRM.Transmittance);
  fprintf('ITM transmission:       %9.4f\n',ifo.Optics.ITM.Transmittance);
  fprintf('PRM transmission:       %9.4f\n',ifo.Optics.PRM.Transmittance);

  figure(10641)
  hndls = loglog(f,sqrt(y1),'-',...         % Quantum Unification  
                 f,sqrt(y10),'-',...        % Seismic
                 f,sqrt(y9),'-',...         % Gravity Gradients
                 f,sqrt(y3),'-',...         % Suspension thermal
                 f,sqrt(y6),'-',...         % Coating Brownian
                 f,sqrt(y11),'--',...       % Coating thermooptic
                 f,sqrt(y5),'--',...        % Substrate brownian
                 f,sqrt(y4),'--',...        % Gas
                 f,sqrt(ys),'k');            % Total Noise
  set(hndls(1:(end)),'LineWidth',5);
  leggravg = strcat('Newtonian background(\beta=',num2str(ifo.Seismic.Beta),')');
  legpower = [num2str(ifo.Laser.Power,'%3.1f') ' W'];
  legend('Quantum noise',...
         'Seismic noise',...
         'Gravity Gradients',...
         'Suspension thermal noise',...
         'Coating Brownian noise',...
         'Coating Thermo-optic noise',...
         'Substrate Brownian noise',...
         'Excess Gas',...
         'Total noise',...
         'Location','NorthEast');
  xlabel('Frequency [Hz]','FontSize',16);
  ylabel('Strain [1/\surdHz]','FontSize',16);
  grid;
  axis([flo fhi 3e-25 1e-21]);
  title(['AdvLIGO Noise Curve: P_{in} = ' legpower],'FontSize',18)  
  
  clrtable=[0.7   0.0   0.9
            0.6   0.4   0.0
            0.0   0.8   0.0
            0.3   0.3   1.0
            1.0   0.2   0.1
            0.0   1.0   0.9
            1.0   0.7   0.0
            0.8   1.0   0.0
            1.0   0.0   0.0
            0.6   0.6   0.6];
  for gag = 1:(length(hndls) - 1)
    set(hndls(gag), 'color',clrtable(gag,:));
  end  
end


if (nargout > 0)
  varargout{1} = 0;
end
switch modeSR
case 0
  nse.ResGas      = y4;
  nse.SuspThermal = y3;
  nse.Quantum     = y1;
  nse.Freq        = f;
  nse.Newtonian   = y9;
  nse.Seismic     = y10;
  nse.Total       = ys;
  nse.MirrorThermal.Total = y2;
  nse.MirrorThermal.SubBrown = y5;
  nse.MirrorThermal.CoatBrown = y6;
  nse.MirrorThermal.SubTE = y8;
  nse.MirrorThermal.CoatTO = y11;
otherwise
  nse.Quantum     = y1;
  nse.Total       = ys;
end
if (nargout > 1)
  varargout{2}    = nse;
end
if (nargout > 2)
  parout.finesse = finesse;
  parout.prfactor = prfactor;
  parout.armpower = armpower;
  parout.bspower = pbs;
  parout.bsthermload = thermalLoad.BS;
  parout.itmthermload = thermalLoad.ITM;
  parout.ifo = ifo;
  
  varargout{3} = parout;
end

% --------------------------------------------------------
% output text

% Report astrophysical scores if so desired
if (makescore == 1)
  sss = int73(nse.Freq, nse.Total, ifo, source);
  sss.Omega = intStoch(nse.Freq, nse.Total, 0, ifo, source);
  if nargout > 0
    varargout{1} = sss;
  end  
end

% Report finesse, power recycling factors
if ( fig > 0 )
  
  disp(sprintf('Finesse:                %7.2f', finesse));
  disp(sprintf('Power Recycling Factor: %7.2f', prfactor))
  disp(sprintf('Arm power:              %7.2f kW', finesse*2/pi*pbs/2/1000));
  disp(sprintf('Power on beam splitter: %7.2f W', pbs))
  PowAbsITM = [finesse*2/pi*ifo.Optics.ITM.CoatingAbsorption/2;...
               ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption/2 ]*pbs;
  M=[ifo.TCS.s_cc,ifo.TCS.s_cs;ifo.TCS.s_cs,ifo.TCS.s_ss];
  S_uncorr=transpose(PowAbsITM)*M*PowAbsITM;
  TCSeff=1-sqrt(ifo.TCS.SRCloss/S_uncorr);
  disp(sprintf('Thermal load on ITM:    %8.3f W', sum(PowAbsITM) ));
  disp(sprintf('Thermal load on BS:     %8.3f W',     ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption    *pbs));
  disp(sprintf(['Reqired TCS efficiency: %8.3f' ...
                '(estimate, see IFOModel.m for definition)'],    TCSeff));  
  if (ifo.Laser.Power*prfactor ~= pbs)
    disp(sprintf('Lensing limited input power: %7.2f W',pbs/prfactor));
  end

  if makescore == 1
     disp(sprintf('BNS Inspiral Range:     %7.2f Mpc',sss.effr0ns))
     disp(sprintf('BBH Inspiral Range:     %7.2f Mpc',sss.effr0bh))
     disp(sprintf('Stochastic Omega:          %4.3g',sss.Omega)) 
  end  
end

return

