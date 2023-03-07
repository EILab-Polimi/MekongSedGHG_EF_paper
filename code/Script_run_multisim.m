%% Script_run_Moltisim

%%% Run a sensitivity analysis on BORG to identify robust dam portfolios for the Mekong that are optimal tradeoffs between
%%% dam sediment trapping, electricity production and GHGs emissions

% addpath('Borg\');

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

par_epsilon = 250;
epsilon =  deltaJ./par_epsilon;  % epsilon for the objectives (1 X n_objectives)
n_obs = length(deltaJ); % number of objectives
NFE = 1E5; % number of function evaluations
n_dec_vars = length(dams);

n_init = 1; % number of initialziations

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

%% Sensitivity Analysis on "GHGs emissions" and "Sediments" values

load('Sediments_Sobol_Analysis');

tic

%run parallel computation
parpool(8);

parfor  iii = 1:n_Sim
    
    dams_SA = dams;

    % % Update GHGS emissions values
    [dams_SA.GHGs_emi_perMWh] = deal(GHGs_emi_perMWh_SA{:,iii});
    [dams_SA.GHGs_emi_perGWh] = deal(GHGs_emi_perGWh_SA{:,iii});
    [dams_SA.GHGs_emi_per_year] = deal(GHGs_emi_per_year_SA{:,iii});
    
    for jjj = 1:size(dams_SA,1)
        
        % % Update Brune Sediment Trapping in each dam (also the one already constucted)
        dams_SA(jjj).SobolMultiplier_brune = SobolSamples(iii, end); %the last sobol parameter indicates the brune efficiency multiplier
          
        % % Update Brune Sediment Trapping in each dam where flushing is possible (varying per country)
        dams_SA(jjj).SobolMultiplier_fl = SobolSamples(iii, (dams_SA(jjj).CountryID + 7)); % the "+7" jumps the first 7 indices that are for the Sediment Yield
        
        if (dams_SA(jjj).Flushing==1 && strcmp(dams_SA(jjj).Status, 'P'))
            
            dams_SA(jjj).TrapEfficiencyBrune = dams_SA(jjj).TrapEfficiencyBrune * dams_SA(jjj).SobolMultiplier_fl * dams_SA(jjj).SobolMultiplier_brune;
            
        elseif dams_SA(jjj).Flushing~=1 %&& strcmp(dams_SA(jjj).Status, 'P'))
            
            dams_SA(jjj).TrapEfficiencyBrune = dams_SA(jjj).TrapEfficiencyBrune * dams_SA(jjj).SobolMultiplier_brune;
            
        end

    end

    % % Update Sediment inputs into each subwatershed
    MS_SA = MS;

    for jjj = 1:size(MS_SA,1)

        if ~isnan(MS_SA(jjj).MekongSubB)
            MS_SA(jjj).SobolMultiplier = SobolSamples(iii, MS(jjj).MekongSubB);
        else
            MS_SA(jjj).SobolMultiplier = 0;
        end

    end

    % % Calculate Sediment inputs into each reach
    QS_in = [MS_SA(:).Yield_tkm2].*[MS_SA(:).directAd].*[MS_SA(:).SobolMultiplier]; % Sediment Yield [t/yr] (local yield [t/km2/yr] * direct area of each node [km2])

    % % Sediment routing
    Theta_S_e_SA = nan(length(QS_in), length(QS_in));

    for jjj = 1:length(QS_in)
        Theta_S_e_SA(jjj, Network.Downstream.Path{jjj,1}{outlet_node}) = QS_in(jjj); % Fractional Sediment flux [kg/yr]
    end
    
    % % Borg
    [decisions, objs, runtime] = borg(n_dec_vars, n_obs, 0, ...
        @(binary_scenario)ObjectiveFunctionMekong_brune_plus_sand(binary_scenario, Network, dams_SA, outlet_node, Theta_S_e_SA), ...
            NFE, epsilon, ExistingDams, upperLimit, parameters);
    
    objs = abs(objs);
    
    summary{iii} = {sparse(round(decisions)); objs; runtime};
    
    parsave_summary('Borg_partial_results_Actual_5000.mat', summary{iii}, iii);
    
end
clear iii

toc

disp('BORG RUN COMPLETED SUCCESSFULLY')

%% Save results

Borg_results.Date = datestr(now,'yy_mm_dd_hh_MM');
Borg_results.epsilon = epsilon;
Borg_results.NFE = NFE;
Borg_results.Nobs = n_obs;
Borg_results.NdecVars = n_dec_vars;
Borg_results.SimulationMode = Simulation_mode;
Borg_results.runtime_min = toc./60;
Borg_results.RESULTS = summary;
Borg_results.Dams = dams;
Borg_results.DamsToExclude = damsToExclude;
Borg_results.ExistingDams = ExistingDams;

save('Borg_results_multisim', 'Borg_results');

disp(['RUN FINISHED AT: ' datestr(now)])
disp(['OUTPUT SAVED AS: ' 'Borg_results_multisim'])
