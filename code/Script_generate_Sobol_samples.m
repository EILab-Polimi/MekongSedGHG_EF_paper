clc; clear;

load('Borg_inputs.mat', 'dams')

n_Cluster = 18;
n_Sim = 1024*5;

% rng
% Sobol_rng_settings = rng;
% save(Sobol_rng_settings);

% load(Sobol_rng_settings);
% rng(Sobol_rng_settings);

Sobol_set = sobolset(n_Cluster);
Sobol_matrix = net(Sobol_set,n_Sim);

%% set Sobol values

minvalues = [repmat(0.2,1,7) repmat(0.1,1,5) repmat(0.5,1,1)];
maxvalues = [repmat(2,1,7) repmat(1,1,5) repmat(1.5,1,1)]; %#ok<REPMAT>

SobolSamples = Sobol_matrix(:,6:end).*(maxvalues-minvalues) + minvalues;

%% dams.GHGs_emi_perMWh

GHGs_emi_perMWh_SA = zeros(length(dams), n_Sim);

for i = 1:length(dams)
    minGHG = dams(i).GHGs_emi_perMWh2_5pct;
    maxGHG = dams(i).GHGs_emi_perMWh97_5pct;
    r = (maxGHG - minGHG).*Sobol_matrix( : , dams(i).Cluster) + minGHG;
    GHGs_emi_perMWh_SA(i,:) = r';
end
clear minGHG; clear maxGHG; clear r; clear i

GHGs_emi_perGWh_SA = GHGs_emi_perMWh_SA * 10^3;

GHGs_emi_per_year_SA = zeros(length(dams), n_Sim);

for i = 1:length(dams)
    GHGs_emi_per_year_SA(i,:) = GHGs_emi_perGWh_SA(i,:) * dams(i).MeanAnnual;
end

%% Save SA matrices

GHGs_emi_perMWh_SA = num2cell(GHGs_emi_perMWh_SA);
GHGs_emi_perGWh_SA = num2cell(GHGs_emi_perGWh_SA);
GHGs_emi_per_year_SA = num2cell(GHGs_emi_per_year_SA);

save('Sediments_Sobol_Analysis.mat','SobolSamples', 'GHGs_emi_perMWh_SA', 'GHGs_emi_perGWh_SA', 'GHGs_emi_per_year_SA');
