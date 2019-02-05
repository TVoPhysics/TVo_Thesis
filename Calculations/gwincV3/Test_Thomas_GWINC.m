source = SourceModel_Thomas;
chirp_mass_old_high_z = [];
horizon_distance_old_high_z = [];
model = IFOModel_Thomas;
chirp_mass_new_high_z = [];
horizon_distance_new_high_z = [];

for i=1:160
%    m = float(i)
    source.NeutronStar.Mass1 = i;
    source.NeutronStar.Mass2 = i;
    ins_mass1 = source.NeutronStar.Mass1;
    ins_mass2 = source.NeutronStar.Mass2;
    tot_mass = ins_mass1 + ins_mass2;
    mu = ins_mass1*ins_mass2/tot_mass;
    M0= mu^(3/5)*tot_mass^(2/5);
    [score, noise] = gwinc_sheila(1,10000,model,source,4);
    horizon_distance_new_high_z(end+1)= score.effr0ns;
    chirp_mass_new_high_z(end+1) = M0;
end

for i=1:160
%    m = float(i)
    source.NeutronStar.Mass1 = i;
    source.NeutronStar.Mass2 = i;
    ins_mass1 = source.NeutronStar.Mass1;
    ins_mass2 = source.NeutronStar.Mass2;
    tot_mass = ins_mass1 + ins_mass2;
    mu = ins_mass1*ins_mass2/tot_mass;
    M0= mu^(3/5)*tot_mass^(2/5);
    [score, noise] = gwinc_Thomas_old(1,10000,model,source,4);
    horizon_distance_old_high_z(end+1)= score.effr0ns;
    chirp_mass_old_high_z(end+1) = M0;
end

  

figure;
loglog(chirp_mass_new_high_z, horizon_distance_new_high_z, 'go-', chirp_mass_old_high_z, horizon_distance_old_high_z, 'bx-');
hold on;
legend('w/ Recursion', 'w/o Recursion');
xlabel('Chirp Mass [Solar Masses]','FontSize',16);
ylabel('Range [MPC]','FontSize',16);
grid on;
