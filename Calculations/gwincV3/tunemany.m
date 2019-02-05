function stat_out = tunemany(varargin)
% TUNEMANY 
%
% Scans over many IFO parameters to find the Gwinc optimum
%
% set doplot=1 to get plots, =2 : more plots

doplot=0;
% Set newdat = 0 to load old data or 1 for new data
newdat = 0;

if nargin == 0
   display('Tuning...')
elseif nargin == 1
  doplot = varargin{1}; % use 1 to make plots, =2 : more plots
elseif nargin == 2
  doplot = varargin{1}; % use 1 to make plots, =2 : more plots
  newdat = varargin{2}; % use 1 for new data
end

f_seis = 9; % Seismic Wall Cutoff frequency
f_Nyq = 8192; % Fake Nyquist frequency

% Range of laser power to search over
Pmin = 3;   % W
Pmax = 125; % W

% Array of laser power to scan over
LaserPower = logspace(log10(Pmin),log10(Pmax),2);

% SRC Detuning phase in radians
Detoon = linspace(0,pi/2,101);

% SRM Transmittance
T_SRM = linspace(0.02,0.3,101);

% ITM Transmittance
T_ITM = linspace(0.005,0.01,2);

noruns = length(LaserPower)*length(Detoon)*length(T_SRM)*length(T_ITM);
nrun = 0;
nave_time = 0;

% Load old data
if newdat == 0
  load tune_many_data

% else run a zillion times and build up the data
else
  for ll = 1:length(T_ITM)
    shotradmode=2;
    for jj = 1:length(LaserPower)
      for kk = 1:length(Detoon)
        for mm = 1:length(T_SRM)

    tic
    [sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,shotradmode,...
                      LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));
    shotradmode=4;
   
    bhs(jj,kk,ll,mm) = sss.effr0bh;
    nhs(jj,kk,ll,mm) = sss.effr0ns;
    omegas(jj,kk,ll,mm) = sss.Omega;

    ntime = toc;
    nrun = nrun + 1;
    nave_time = nave_time + ntime;
    ave_run_time = nave_time / nrun;
    estimated_time = (noruns - nrun) * ave_run_time;


        end
      end
     disp([num2str(estimated_time/60) ' minutes remaining'])
    end
    save tune_many_data LaserPower Detoon T_ITM T_SRM bhs nhs omegas
  end  

  display(' ')
  display('Done generating the gwinc hypercube')
  display(' ')

end



if doplot > 0

  % New smoov hot map
  myhot = hot(3000);


% ------------- Black Holes -----------------------------------------------
figure(11)
[maxx,imaxx] = max(bhs(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(bhs),imaxx); % get multi-D indices from flat index

fprintf('\n)ptimized for BH-BH (30/30)\n');
fprintf('-----------------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for 30/30 Inspirals (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
grid minor
print -dpng BHTune2.png
save noise_BHBH2 nnn

try
figure(12)
surf(LaserPower,Detoon*180/pi,bhs(:,:,ll,mm)')
view(0,90)
axis tight
%axis([1 max(LaserPower) 0 max(Detoon*180/pi) 0 max(max(bhs))])
set(gca,'XScale','log')
shading interp
colormap(myhot)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('BH-BH (30/30) Inspiral Range [Mpc]')
colorbar
print -dpng BHsurf2.png
end

if doplot > 1
  figure(13)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(14)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end

% ----------------- Neutron Stars ----------------------------------------
figure(21)
[maxx,imaxx] = max(nhs(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(nhs),imaxx); % get multi-D indices from flat index

fprintf('\nOoptimized for NS-NS\n');
fprintf('---------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for NS/NS Inspirals (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
print -dpng NSTune2.png
save noise_NSNS2 nnn

try
figure(22)
surf(LaserPower,Detoon*180/pi,nhs(:,:,ll,mm)')
view(0,90)
axis tight
set(gca,'XScale','log')
shading interp
colormap(myhot)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('NS-NS Inspiral Range [Mpc]')
colorbar
print -dpng NSsurf2.png
end

if doplot > 1
  figure(23)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(24)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end

% -------------- Big Bang Leftovers --------------------------------------
figure(31)
[maxx,imaxx] = min(omegas(:));    % Get max BHS score for flattened bhs variable
[jj,kk,ll,mm] = ind2sub(size(omegas),imaxx); % get multi-D indices from flat index

fprintf('\nOptimized for Stochastic\n');
fprintf('--------------------------------\n');

[sss,nnn] = gwinc(f_seis,f_Nyq,IFOModel,SourceModel,3,...
               LaserPower(jj),Detoon(kk),T_SRM(mm),T_ITM(ll));

title(['Optimized for Stochastic (P= ' num2str(LaserPower(jj),3) ...
       ' W, \phi_{SRM} = ' num2str(Detoon(kk)*180/pi,2) ' deg)'])
print -dpng StochTune2.png

try
figure(32)
surf(LaserPower,Detoon*180/pi,-log10(omegas(:,:,ll,mm)'))
view(0,90)
axis tight
set(gca,'XScale','log')
shading interp
colormap(myhot)
brighten(0.3)
xlabel('Laser Power [W]')
ylabel('SRM Detune [deg]')
title('-log_{10}[Omega_{GW}]')
colorbar
print -dpng Stochsurf2.png
end

if doplot > 1
  figure(33)
  surf(Detoon*180/pi,T_SRM,reshape(nhs(jj,:,ll,:),length(Detoon),length(T_SRM))')
  view(0,90)
  axis tight
  shading interp
  colormap(myhot)
  xlabel('SRM Detune [deg]')
  ylabel('SRM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar

  figure(34)
  surf(LaserPower,T_ITM,reshape(nhs(:,kk,:,mm),length(LaserPower),length(T_ITM))')
  view(0,90)
  axis tight
  set(gca,'XScale','log')
  shading interp
  colormap(myhot)
  xlabel('Laser Power [W]')
  ylabel('ITM Transmission')
  title('NS-NS Inspiral Range [Mpc]')
  colorbar
end



end

stat_out = newdat;


