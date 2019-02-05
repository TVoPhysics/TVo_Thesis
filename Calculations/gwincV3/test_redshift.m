source = SourceModel;
ifo = IFOModel;
f = logspace(log10(1),log10(10000),3000);
[score,noise]=gwinc(1,10000,ifo,source,2);
h2 = (noise.Total)*.1;

chirp_mass_orig = [];
horizon_distance_orig = [];

chirp_mass_sheila = [];
horizon_distance_sheila = [];

chirp_mass_thomas = [];
horizon_distance_thomas = [];

source.Distance = 1500

tic
for i=1:100 
    source.NeutronStar.Mass1 = i;
    source.NeutronStar.Mass2 = i;
    ins_mass1 = source.NeutronStar.Mass1;
    ins_mass2 = source.NeutronStar.Mass2;
    tot_mass = ins_mass1 + ins_mass2;
    mu = ins_mass1*ins_mass2/tot_mass;
    M0= mu^(3/5)*tot_mass^(2/5);
    score = int73_Thomas_new1(f, h2, ifo, source);
    horizon_distance_thomas(end+1)= score.effr0ns; 
    chirp_mass_thomas(end+1) = score.m_ns;
end
toc

tic
for i=1:100
    source.NeutronStar.Mass1 = i;
    source.NeutronStar.Mass2 = i;
    ins_mass1 = source.NeutronStar.Mass1;
    ins_mass2 = source.NeutronStar.Mass2;
    tot_mass = ins_mass1 + ins_mass2;
    mu = ins_mass1*ins_mass2/tot_mass;
    M0= mu^(3/5)*tot_mass^(2/5);
    score = int73(f, h2, ifo, source);
    horizon_distance_orig(end+1)= score.effr0ns*2.2643;
    chirp_mass_orig(end+1) = M0;
end
toc

for i=1:100
    source.NeutronStar.Mass1 = i;
    source.NeutronStar.Mass2 = i;
    ins_mass1 = source.NeutronStar.Mass1;
    ins_mass2 = source.NeutronStar.Mass2;
    tot_mass = ins_mass1 + ins_mass2;
    mu = ins_mass1*ins_mass2/tot_mass;
    M0= mu^(3/5)*tot_mass^(2/5);
    score = int73_noredshift(f, h2, ifo, source.NeutronStar.Mass1);
    horizon_distance_sheila(end+1)= score.effr0;
    chirp_mass_sheila(end+1) = M0;
end

dist_z = lumDist2Z(horizon_distance_sheila/1000);
chirp_mass_sheila_new = chirp_mass_sheila ./ (1+dist_z);

figure;
loglog(chirp_mass_sheila_new, horizon_distance_sheila, 'go-', ....
chirp_mass_orig, horizon_distance_orig, 'bx-', chirp_mass_thomas, ... 
horizon_distance_thomas, 'rx-');
%loglog(chirp_mass_thomas,horizon_distance_thomas, 'rx-');
hold on;
legend('Sheila', 'Original', 'Thomas');
xlabel('Chirp Mass [Solar Masses]','FontSize',16);
ylabel('Horizon Distance [MPC]','FontSize',16);
grid on;