function score = int73(f, h2, ifo, source)
% score = int73(f, h2, ifo, source)
% Takes the IFO noise and the frequency grid and make the inspiral
% score via the f^7/3 integral
% includes a conversion to horizon distance with a factor of 2.2643 in the
% getOverlapforMassPair function
% also includes a recusive algorithm to compute the correct chirpmass to horizon distance 
% relation that includes redshift.
 
  % try to find global f_bounce, set by gwinc in suspR
  % this is used to define the start of the integral
  global f_bounce
  if isempty(f_bounce)
    f_bounce = 9.1416;
  end

  % Constants
  MSol = ifo.Constants.MSol;
  Mpc = ifo.Constants.Mpc;
  c = ifo.Constants.c;
  g = pi*MSol/c;


  % Start the integration 0.1 Hz above the bounce mode
  n = find(f > (f_bounce+0.1));
  f = f(n);
  h2 = h2(n);

  integrand73 = (g^2) ./ ((g*f).^(7/3) .* h2);
  x73all = cumtrapz(f,integrand73);

  x73 = x73all(end);
  score.r0 = sqrt(x73*(3/20)^(5/3)*(5/(192*pi))); % Finn 96 def.
  % volume to distance conversion; from P Sutton; see T1200351
  score.r0 = (3*1.8375)^(1/3)*score.r0*MSol/Mpc;   % adjust for survey volume (in Mpc)

  
  % Build arrays to do quick distance to redshift conversion
  d_list = [.1:.1:100];
  red_list = lumDist2Z(d_list);
     
  accuracy = .001; %the iteration will work to within a ~ 2.5 MPC
  % 1.4/1.4 Binary Neutron Stars -------------------------------------------------
  % Initial estimates from reference source file
  ins_mass1_ns = source.NeutronStar.Mass1*MSol;
  ins_mass2_ns = source.NeutronStar.Mass2*MSol;
  tot_mass_ns = ins_mass1_ns + ins_mass2_ns;
  mu_ns = ins_mass1_ns*ins_mass2_ns/tot_mass_ns; %reference reduced mass
  M0_ns= mu_ns^(3/5)*tot_mass_ns^(2/5); %reference chirp mass   
  distance_ns = source.NeutronStar.Distance;   % in Mc 
  
  % Recursion to get an accurate effective range and red-shifted chirp mass
  m_obs_ns = M0_ns; 
  tot_mass_obs_ns=tot_mass_ns;
  z_test = 0;
  r0 = score.r0;
  
  for i = 1:500;
%  if false
%      [effr0ns, z, m_obs_ns,tot_mass_obs_ns] = CalcRedshift(m_obs_ns, M0_ns, tot_mass_ns, distance_ns, f, x73all, MSol, c, d_list, red_list,r0);  
%  else
      epsilon=MSol*5e-3;
      mold=m_obs_ns;
      [effr0ns, z, m_obs_ns,tot_mass_obs_ns] = CalcRedshift(m_obs_ns, M0_ns, tot_mass_ns, distance_ns, f, x73all, MSol, c, d_list, red_list,r0);  
      [effr0ns2, z2, m_obs_ns2,tot_mass_obs_ns2] = CalcRedshift(m_obs_ns+epsilon, M0_ns, tot_mass_ns, distance_ns, f, x73all, MSol, c, d_list, red_list,r0);
      ff=m_obs_ns-mold; 
      ff2=m_obs_ns2-mold-epsilon;
      m_obs_ns=m_obs_ns-(ff)/(ff2-ff)*epsilon;
%  end
  
  error = abs(z-z_test);
  if error < accuracy  || isnan(error)==1; % test how close the function is to converting
      break
  else
       z_test = z; 
       continue
  end
  end
        disp([num2str(i),' ',num2str(M0_ns/MSol)]);
  score.effr0ns = effr0ns;
  score.m_ns = M0_ns/MSol;
%  score.ff = ff;

 %{
  % Black Holes  -------------------------------------------------
  ins_mass1_bh = source.BlackHole.Mass1*MSol;
  ins_mass2_bh = source.BlackHole.Mass2*MSol;
  tot_mass_bh = ins_mass1_bh + ins_mass2_bh;
  mu_bh = ins_mass1_bh*ins_mass2_bh/tot_mass_bh;
  M0_bh= mu_bh^(3/5)*tot_mass_bh^(2/5);   
  distance_bh = source.BlackHole.Distance;   % in Mpc
  
  % Recursion to get an accurate effective range and red-shifted chirp mass
  m_obs_bh = M0_bh; 
  tot_mass_obs_bh=tot_mass_bh;
  
  effr0bh_test = 0;
 for i = 1:10;
  [r_bh,m_bh] = getOverlapForMassPair(m_obs_bh, tot_mass_obs_bh, distance_bh, f, x73all, MSol, c, d_list, red_list);
  score.effr0bh = r0 * r_bh; % Output: horizon distance in MPC
  z = interp1(d_list,red_list,score.effr0bh/1000,'nearest');
  m_obs_bh = M0_bh * (1+z);
  tot_mass_obs = tot_mass_bh * (1+z);
  score.m_bh = M0_bh/MSol; %Redshift of Mass
  error = abs(score.effr0bh-effr0bh_test); 
    if error < accuracy  || isnan(error)==1;
      break
  else
      continue
    end
  effr0bh_test = score.effr0bh;
 end
%}
end

%  ------------------------------------------------- --------------------------
%
% return the score in Mpc for a inspiral source
% the source struct should contain
%   Mass1 and Mass2 in units?
%   Distance in units?
%
% x73all is the cumulative range integral

function [effr0, z, m_obs,tot_mass_obs] = CalcRedshift(m_obs, M0,tot_mass_ref, distance, f, x73all, MSol, c, d_list, red_list,r0)

  tot_mass_obs = m_obs/M0*tot_mass_ref;
  [r,m] = getOverlapForMassPair(m_obs, tot_mass_obs, distance, f, x73all, MSol, c, d_list, red_list);
  effr0 = r0 * r; % Output: horizon distance in MPC, score.r0 are the constants in line 35
  z = interp1(d_list,red_list,effr0/1000); % Use luminosity to redshift conversion
  m_obs = M0 * (1+z); % Update the observed chirpmass by a redshift factor, this will be re-inputted into the next iteration

end

function [score,M0] = getOverlapForMassPair(M0, tot_mass, distance, f, x73all, MSol, c, d_list, red_list)
    
  z = 0;
  MM = M0;
    
  f2max = 2 * (1980/(z+1)) * (MSol/tot_mass);
% if (f2max < f(1))
%    error('f_isco is less than f_low')
% end

  x73 = x73all(end);
  zetafmax = interp1(f, x73all, f2max,'nearest') / x73;

  % debug: make sure point lies on curve
  %semilogx(f, x73all, f2max, x73 * zetafmax, 'x');
  %pause
  
  score = sqrt(zetafmax) * (MM/1.2188/MSol)^(5/6) * 2.2643; % 2.2643 is conversion to horizon distance

end
