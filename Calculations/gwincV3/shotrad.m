function n = shotrad(f, ifo)
% Quantum noise model
% All references to Buonanno & Chen PRD 64 042006 (2001) (hereafter BnC)
% Updated to include losses DEC 2006 Kirk McKenzie using BnC notation

% f                                            % Signal Freq. [Hz]
lambda  = ifo.Laser.Wavelength;               % Laser Wavelength [m]
hbar    = ifo.Constants.hbar;                 % Plancks Constant [Js]
c       = ifo.Constants.c;                    % SOL [m/s]
Omega   = 2*pi*f;                             % [BnC, table 1] Signal angular frequency [rads/s]
omega_0 = 2*pi*c/lambda;                      % [BnC, table 1] Laser angular frequency [rads/s]

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
L       = ifo.Infrastructure.Length;          % Length of arm cavities [m]
l       = ifo.Optics.SRM.CavityLength;        % SRC Length [m]
T       = ifo.Optics.ITM.Transmittance;       % ITM Transmittance [Power]
m       = ifo.Materials.MirrorMass;            % Mirror mass [kg]

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bsloss  = ifo.Optics.BSLoss;                  % BS Loss [Power]
mismatch = 1 - ifo.Optics.coupling;           % Mismatch
mismatch = mismatch + ifo.TCS.SRCloss;        % Mismatch

% BSloss + mismatch has been incorporated into a SRC Loss
lambda_SR = mismatch + bsloss;                % SR cavity loss [Power]
lambda_PD = 1 - ifo.Optics.PhotoDetectorEfficiency; % Loss in PD Process [Power]

tau     = sqrt(ifo.Optics.SRM.Transmittance); % SRM Transmittance [amplitude]
rho     = sqrt(1 - tau^2 - lambda_SR);        % SRM Reflectivity [amplitude]

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds      = ifo.Optics.SRM.Tunephase;           % SRC Detunning
eta     = ifo.Optics.Quadrature.dc;           % Homodyne Readout phase

phi     = (pi-ds)/2;                          % [BnC, between 2.14 & 2.15] SR Detuning
Phi     = mod(l*Omega/c,2*pi);                % [BnC, between 2.14 & 2.15] SRC one pass phase shift

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lambda_arm = ifo.Optics.Loss*2;               % [BnC, after 5.2] Round Trip loss in arm [Power]
gamma_ac = T*c/(4*L);                         % [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
epsilon = lambda_arm/(2*gamma_ac*L/c);        % [BnC, after 5.2] Loss coefficent for arm cavity
Epsilon = 2*epsilon./(1+(Omega/gamma_ac).^2); % [BnC, 5.3] Loss Parameter

I_0     = ifo.gwinc.pbs;                      % [BnC, Table 1] Power at BS (Power*prfactor) [W]
I_SQL   = (m*L^2*gamma_ac^4)/(4*omega_0);     % [BnC, 2.14] Power to reach free mass SQL
Kappa  = 2*((I_0/I_SQL)*gamma_ac^4)./...
         (Omega.^2.*(gamma_ac^2+Omega.^2));   % [BnC 2.13] Effective Radiation Pressure Coupling
beta    = atan(Omega./gamma_ac);              % [BnC, after 2.11] Phase shift of GW SB in arm
h_SQL   = sqrt(8*hbar./(m*(Omega*L).^2));     % [BnC, 2.12] SQL Strain

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Coefficients [BnC, Equations 5.8 to 5.12]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C11_L   = sqrt(1-lambda_PD) * ( (1+rho^2) * ( cos(2*phi) + Kappa/2 * sin(2*phi) ) -...
          2*rho*cos(2*beta) - 1/4*epsilon * ( -2 * (1+exp(2*i*beta)).^2 * rho + 4 * (1+rho^2) *...
          cos(beta).^2*cos(2*phi) + ( 3+exp(i*2*beta) ) .* Kappa * (1+rho.^2) * sin(2*phi) ) + ...
        lambda_SR * ( exp(2*i*beta)*rho-1/2 * (1+rho^2) * ( cos(2*phi)+Kappa/2 * sin(2*phi) ) ) );

C22_L   = C11_L;

C12_L   = sqrt(1-lambda_PD) * tau^2 * ( - ( sin(2*phi) + Kappa*sin(phi).^2 )+...
 1/2*epsilon*sin(phi) * ( (3+exp(2*i*beta)) .* Kappa * sin(phi) + 4*cos(beta) .^2 * cos(phi)) + ...
 1/2*lambda_SR * ( sin(2*phi)+Kappa*sin(phi).^2) );

C21_L   = sqrt(1-lambda_PD) * tau^2 * ( (sin(2*phi)-Kappa*cos(phi).^2 ) + ...
     1/2*epsilon*cos(phi) * ( (3+exp(2*i*beta) ).*Kappa*sin(phi) - 4*cos(beta).^2*sin(phi) ) + ...
     1/2*lambda_SR * ( -sin(2*phi) + Kappa*cos(phi).^2) );

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
D1_L    = sqrt(1-lambda_PD) * ( - (1+rho*exp(2*i*beta) ) * sin(phi) + ...
        1/4*epsilon * ( 3+rho+2*rho*exp(4*i*beta) + exp(2*i*beta)*(1+5*rho) ) * sin(phi)+ ...
        1/2*lambda_SR * exp(2*i*beta) * rho * sin(phi) );

D2_L    = sqrt(1-lambda_PD) * ( - (-1+rho*exp(2*i*beta) ) * cos(phi) + ...
        1/4*epsilon * ( -3+rho+2*rho*exp(4*i*beta) + exp(2*i*beta) * (-1+5*rho) ) * cos(phi)+ ...
        1/2*lambda_SR * exp(2*i*beta) * rho * cos(phi) );

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
P11     = 1/2*sqrt(1-lambda_PD) * sqrt(lambda_SR) * tau *...
          ( -2*rho*exp(2*i*beta)+2*cos(2*phi)+Kappa*sin(2*phi) );
P22     = P11;
P12     = -sqrt(1-lambda_PD) * sqrt(lambda_SR)*tau*sin(phi)*(2*cos(phi)+Kappa*sin(phi) );
P21     =  sqrt(1-lambda_PD) * sqrt(lambda_SR)*tau*cos(phi)*(2*sin(phi)-Kappa*cos(phi) );

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Q11     = sqrt(lambda_PD) *...
          ( exp(-2*i*beta)+rho^2*exp(2*i*beta)-rho*(2*cos(2*phi)+Kappa*sin(2*phi)) + ...
          1/2*epsilon*rho * (exp(-2*i*beta)*cos(2*phi)+exp(2*i*beta).*...
                            ( -2*rho-2*rho*cos(2*beta)+cos(2*phi)+Kappa*sin(2*phi) ) + ...
        2*cos(2*phi)+3*Kappa*sin(2*phi))-1/2*lambda_SR*rho *...
            ( 2*rho*exp(2*i*beta)-2*cos(2*phi)-Kappa*sin(2*phi) ) );
Q22     = Q11;
Q12     = 0;
Q21     = 0;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N11     = sqrt(1-lambda_PD) * sqrt(epsilon/2)*tau *(Kappa.*(1+rho*exp(2*i*beta))*sin(phi)+...
          2*cos(beta).*(exp(-i*beta)*cos(phi)-rho*exp(i*beta).*(cos(phi)+Kappa*sin(phi))));
N22     = -sqrt(1-lambda_PD)*sqrt(2*epsilon)*tau*(-exp(-i*beta)+rho*exp(i*beta)).*...
           cos(beta)*cos(phi);
N12     = -sqrt(1-lambda_PD)*sqrt(2*epsilon)*tau*(exp(-i*beta)+rho*exp(i*beta)).*...
           cos(beta)*sin(phi);
N21     = sqrt(1-lambda_PD)*sqrt(2*epsilon)*tau*(-Kappa*(1+rho)*cos(phi)+...
          2*cos(beta).*(exp(-i*beta)+rho*exp(i*beta)).*cos(beta)*sin(phi));


%>>>>>>>>    QUANTUM NOISE POWER SPECTRAL DENSITY [BnC, 5.13]   <<<<<<<<<<<<<<<<<

n = h_SQL.^2./(2*Kappa*tau^2.*abs(D1_L*sin(eta)+D2_L*cos(eta)).^2).*(...
    abs(C11_L*sin(eta)+C21_L*cos(eta)).^2+...
    abs(C12_L*sin(eta)+C22_L*cos(eta)).^2+abs(P11*sin(eta)+P21*cos(eta)).^2+...
    abs(P12*sin(eta)+P22*cos(eta)).^2+abs(Q11*sin(eta)+Q21*cos(eta)).^2+...
    abs(Q12*sin(eta)+Q22*cos(eta)).^2+abs(N11*sin(eta)+N21*cos(eta)).^2+...
    abs(N12*sin(eta)+N22*cos(eta)).^2);

