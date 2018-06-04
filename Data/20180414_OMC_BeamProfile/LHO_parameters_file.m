
% LHO parameters file
ISI2OM1 = 1.59; % edge of ISI table next to HAM5 viewport to OM1
SRM2OM1 = 3.4206; % from Zemax
OM1S = 1/4.6; % T1200410
OM1toOM2 = 1.39; % matches what we measured
OM2S = 1/(1.7); % T1200410
OM2toWaist = 0.64+0.456;
%SRMtoWaist = M5*M4*M3*M2*M1;
%WaisttoSRM = M1*M2*M3*M4*M5;

% Transmission through SRM
SRMROC = -5.678;
SRMthick = 74.73E-3;
n1 = 1; n2 = 1.4496;

% through SRC
SRMtoSR2 = 15.739;
SR2ROC = -6.425; 
SR2toSR3 = 15.461;
SR3ROC = 36.013; % galaxy
SR3toBS = 24.368;
BStoITM = 0;
ITMROC = -1939;
