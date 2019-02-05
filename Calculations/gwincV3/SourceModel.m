function source = SourceModel(varargin);

% SOURCEMODEL returns a structure describing gravity wave sources for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 
% 


%% Neutron Star
source.NeutronStar.Mass1 = 1.4;	% Solar Mass
source.NeutronStar.Mass2 = 1.4; % Solar Mass
source.NeutronStar.Distance = 23; % Megaparces

%% Black Hole
source.BlackHole.Mass1 = 30; % Solar Mass
source.BlackHole.Mass2 = 30; % Solar Mass
source.BlackHole.Distance = 100; % Megaparces

%% NS/BH 1.4/30 Short Hard Gamma Ray Bursts
source.SHGRB.Mass1 = 1.4; % Solar Mass
source.SHGRB.Mass2 = 30; % Solar Mass
source.SHGRB.Distance = 40; % Megaparces

%% Stochastic
source.Stochastic.powerlaw = 0;           %
source.Stochastic.integration_time = 1; % years
source.Stochastic.Hubble = 72;          % Hubble constant in km/s/Mpc
source.Stochastic.confidence = 0.9;    % confidence level for frequentist upper limit

%% Pulsar

source.Pulsar.name1 = 'Crab Pulsar';      % First pulsar
source.Pulsar.distance1 = 1.9;            % kiloparsecs          % kiloparsecs
source.Pulsar.I31 = 3e38;                 % Moment of inertia, in mks units
source.Pulsar.rotation_frequency1 = 29.8; % rotational frequency
source.Pulsar.frequency_multiple1 = 2;    % frequency multiplier for GW emission
source.Pulsar.integration_time1 = 1;      % integraton time in years

source.Pulsar.name2 = 'Sco X-1';          % Second pulsar
source.Pulsar.distance2 = 2.8;            % kiloparsecs
source.Pulsar.I32 = 1e38;                 % Moment of inertia, in mks units
source.Pulsar.rotation_frequency2 = 310;  % rotational frequency
source.Pulsar.frequency_multiple2 = 2;    % frequency multiplier for GW emission
source.Pulsar.integration_time2 = 1;      % integraton time in years

%%%%%%

return
