
addpath(genpath(pwd))

load('Borg_inputs');

%% Set the simulation mode and decide if to exclude any dam

Simulation_mode = 'Actual'; % use 'Pristine', 'Actual'
% % Pristine: Situation with no dams
% % Actual: Current situation (dams with status "built" or "under construction" cannot be excluded)

damsToExclude = zeros(length(dams),1); % damsToExclude is 1 for each dam that should not be included in the portfolio analysis

if strcmp(Simulation_mode, 'Pristine')
    ExistingDams = zeros(length(dams),1);
elseif strcmp(Simulation_mode, 'Actual')
    ExistingDams = zeros(length(dams),1);
    ExistingDams(ismember({dams.Status}, {'C','E'})) = 1;
end

%% Set Borg parameters

binary_scenario = ones(length(dams), 1);
JAllDams = ObjectiveFunctionMekong_brune_plus_sand(binary_scenario, Network, dams, outlet_node, Theta_S_e); % objective with all dams

binary_scenario = zeros(length(dams), 1);
JNoDams = ObjectiveFunctionMekong_brune_plus_sand(binary_scenario, Network, dams, outlet_node, Theta_S_e);  % objective with no dams

deltaJ = abs(JAllDams - JNoDams);
epsilon =  deltaJ./250;  % epsilon for the objectives (1 X n_objectives)
n_obs = length(deltaJ); % number of objectives
NFE = 1E5; % number of function evaluations
n_dec_vars = length(dams);

n_init = 1; % number of initialiations

randseed = 28;
parameters = {'rngstate', randseed}; % seed

%% Inlude existing dams and/or exclude certain dam sites
% % However, if the lower and upper limit are the same in Borg,
% % the algorithm seems to set the decision variable to nan, causing the objective function to crash
% % This is avoided by setting dams that are build before the date cutoff to 0.9 which will then be rounded to 1 by the algorithm

damsToExclude(ExistingDams == 1) = 0; % dams that already exist cannot be excluded

%ExistingDams is also the lowerLimit
ExistingDams(ExistingDams == 1) = 0.9;

% % Upper limit for decision variables
upperLimit = ones(n_dec_vars,1);
if exist('damsToExclude')
    upperLimit(find(damsToExclude)) = 0.1; % set the upper limit of the decision variable to 0.1 (i.e., 0) for all dams that should be excluded
end

%% RUN Borg
tic

% % Calculate Sediment inputs into each reach
QS_in = [MS(:).Yield_tkm2].*[MS(:).directAd]; % Sediment Yield [t/yr] (local yield [t/km2/yr] * direct area of each node [km2])

% % Sediment routing
Theta_S_e_SA = nan(length(QS_in), length(QS_in));

for jjj = 1:length(QS_in)
    Theta_S_e_SA(jjj, Network.Downstream.Path{jjj,1}{outlet_node}) = QS_in(jjj); % Fractional Sediment flux [kg/yr]
end

% % Borg
[decisions, objs, runtime] = borg(n_dec_vars, n_obs, 0, ...
    @(binary_scenario)ObjectiveFunctionMekong_brune_plus_sand(binary_scenario, Network, dams, outlet_node, Theta_S_e_SA), ...
        NFE, epsilon, ExistingDams, upperLimit, parameters);

objs = abs(objs);

results_avgscenario = {round(decisions); objs; runtime};

toc

disp('BORG RUN COMPLETED SUCCESSFULLY')

%% save results

Borg_results_avg.Date = datestr(now,'yy_mm_dd_hh_MM');
Borg_results_avg.epsilon = epsilon;
Borg_results_avg.NFE = NFE;
Borg_results_avg.Nobs = n_obs;
Borg_results_avg.NdecVars = n_dec_vars;
Borg_results_avg.SimulationMode = Simulation_mode;
Borg_results_avg.runtime_min = toc./60;
Borg_results_avg.RESULTS = results_avgscenario;
Borg_results_avg.Dams = dams;
Borg_results_avg.DamsToExclude = damsToExclude;
Borg_results_avg.ExistingDams = ExistingDams;

save('Borg_results_averagesim', 'Borg_results_avg');
