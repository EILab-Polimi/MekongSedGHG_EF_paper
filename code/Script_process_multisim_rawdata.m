%% load data 

load('Borg_inputs');

load('Borg_results_multisim_raw.mat', 'Borg_results')

%% extract planned dams and calculate distance from outlet

% "Actual" simulation mode -> only consider "Planned" dams
outlet = 128;
index_dams_P = [];
index_dams_E = [];

for iii = 1:length(dams)
    
    dams(iii).dist2out = Network.Downstream.Distance{dams(iii).FromNode}(outlet);
    if strcmp(dams(iii).Status, 'P')
        index_dams_P = [index_dams_P; iii];
    else
        index_dams_E = [index_dams_E; iii];
    end
end
clear iii

dams_included = dams(index_dams_P);

%% calculate total average probability of inclusion

% n simulation sensitivity analysis
n_Sim = length(Borg_results.RESULTS);

%define energy generation scenario via a range
range_HPP = 1.3e5:10e3:2.7e5;

% define matrix with probability of inclusion for each dam and each energy
% generation scenario

prob_dam_inclusion = zeros(length(index_dams_P),length(range_HPP)-1);
dam_pos = zeros(size(index_dams_P));

% define particular percentage to be found

pr_detail = [0.1 0.5 0.9];

% HP_scenario_percentage_pos contains for each considered dam (row) and percentage in pr_detail (column) the position in range where the dam falls more closely 
HP_scenario_percentage = zeros(length(index_dams_P),length(pr_detail));

wb = waitbar(1/length(index_dams_P), ['dam 1/' num2str(length(index_dams_P))]); %open waitbar

for d = 1:length(index_dams_P)
    
    waitbar(d/length(index_dams_P), wb, ['dam ' num2str(d) '/' num2str(length(index_dams_P)) ]); % update waitbar

    dam_pos = index_dams_P(d);
    range_change = range_HPP;
    
    for r = 1:length(range_HPP)-1
        
        range_pos = [range_HPP(r) range_HPP(r+1)];

        count_dams = 0; %count the total number of times the dam d appears in the PO solutions within range r 
        count_POsolutioninrange(r) = 0; %count the total number of PO solutions with HPproduction within range r
        
        for i = 1:n_Sim
            
            Portfolios = full(Borg_results.RESULTS{1, i}{1, 1});
            
            JHP_i = Borg_results.RESULTS{1, i}{2, 1}(:,2);
            n_POsolutioninrange = and(JHP_i>range_pos(1), JHP_i<range_pos(2)); %number of PO solution with HPproduction within range r
            
            count_POsolutioninrange(r) = count_POsolutioninrange(r) + sum(n_POsolutioninrange);
            
            count_dams = count_dams + sum(Portfolios(n_POsolutioninrange,dam_pos));
            
        end
        
        prob_dam_inclusion(d,r) = count_dams/ count_POsolutioninrange(r);
        
    end
    
    range_change(isnan(prob_dam_inclusion(d,:))) = [];
    
    for s=1:length(pr_detail)
        
        ix = find(prob_dam_inclusion(d,:)>pr_detail(s),1);
            
        if ix ==1
            HP_scenario_percentage(d,s) = range_change(ix);
        elseif isempty(ix)
            HP_scenario_percentage(d,s) = range_change(end);
        else
            HP_scenario_percentage(d,s) = interp1(prob_dam_inclusion(d,ix-1:ix),range_change(ix-1:ix), pr_detail(s) );
        end                
    end
end

close(wb)

%% calculate probability of inclusion for each sensitivity scenario

%chose dams
OBJECTID = [dams.OBJECTID]';

% n simulation sensitivity analysis
n_Sim = length(Borg_results.RESULTS);

% define particular percentage to be found

pr_detail = [0.1 0.5 0.9];
HP_scenario_percentage_allsim = cell(length(index_dams_P),1);

wb = waitbar(1/length(index_dams_P), ['dam 1/' num2str(length(index_dams_P))]); %open waitbar

for d = 1:length(index_dams_P)
    
    waitbar(d/length(index_dams_P), wb, ['dam ' num2str(d) '/' num2str(length(index_dams_P)) ]); % update waitbar
    
    [HP_scenario_percentage_allsim{d,:}] = deal(zeros(length(pr_detail),n_Sim));
    
    dam_pos = index_dams_P(d);
    
    for i = 1:n_Sim %sensitivity scenario 
        
        n_POsolutions = size(Borg_results.RESULTS{1, i}{2, 1}(:,2),1); %n PO portfolios in n_Sim
            
        prob_dam_inclusion_sim = zeros(1,length(range_HPP)-1);

        for r = 1:length(range_HPP)-1

            range_change = range_HPP;
            
            range_pos = [range_HPP(r) range_HPP(r+1)];

            Portfolios = full(Borg_results.RESULTS{1, i}{1, 1});
            
            JHP_i = Borg_results.RESULTS{1, i}{2, 1}(:,2);
            n_POsolutioninrange = and(JHP_i>range_pos(1), JHP_i<range_pos(2)); %n of PO solutions with HPP in range r
            
            count_POsolutioninrange = sum(n_POsolutioninrange); %count the total number of times the dam d appears within the range r   
            
            count_dams = sum(Portfolios(n_POsolutioninrange,dam_pos));
            
            prob_dam_inclusion_sim(r) = count_dams/count_POsolutioninrange;
            
            range_change(isnan(prob_dam_inclusion_sim)) = [];
            prob_dam_inclusion_sim(isnan(prob_dam_inclusion_sim)) = [];
            
        end
        
        for s=1:length(pr_detail)
            
            ix = find(prob_dam_inclusion_sim>pr_detail(s),1);
            
            if ix ==1
                HP_scenario_percentage_allsim{d}(s,i) = range_change(ix);
            elseif isempty(ix)
                HP_scenario_percentage_allsim{d}(s,i) = range_change(end);
            else
                                         
                HP_scenario_percentage_allsim{d}(s,i) = interp1(prob_dam_inclusion_sim(ix-1:ix),range_change(ix-1:ix), pr_detail(s) );
            end
        end
    
    end
    
end

close(wb)

%% save results

save('Borg_results_multisim' , 'HP_scenario_percentage', 'HP_scenario_percentage_allsim', 'pr_detail', 'prob_dam_inclusion', 'range_HPP');