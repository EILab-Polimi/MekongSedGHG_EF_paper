%% Figure 2

% load data
addpath(genpath(pwd))

load('Borg_results_averagesim.mat')
load('Borg_inputs');

% "Actual" simulation mode -> only consider "Planned" dams
outlet_node = 128;
index_dams_P_bin = zeros(size(dams))';
index_dams_E_bin = zeros(size(dams))';

for iii = 1:length(dams)
    
    dams(iii).dist2out = Network.Downstream.Distance{dams(iii).FromNode}(outlet_node);
    if strcmp(dams(iii).Status, 'P')
        index_dams_P_bin(iii) = 1;
    else
        index_dams_E_bin(iii) = 1;
    end
end
clear iii

index_dams_P_bin = index_dams_P_bin ==1;
index_dams_E_bin = index_dams_E_bin ==1;


dams_included = dams(index_dams_P_bin);

Portfolios = Borg_results_avg.RESULTS{1, 1}; % displays the matrix in the standard way (with 0)
J = Borg_results_avg.RESULTS{2, 1};

[J, ord] = sortrows(J,2); % sorts the rows of obj_values according to increasing value of J2(HPP) and create the relative order of sorting variable ("ord")
J1_sed = J(:,1)/10^6; J2_HPP = J(:,2); J3_GHGs = J(:,3);
Portfolios = Portfolios(ord,:); % sorts the portfolios according to increasing value of the correspondent J2(HPP)

% Figure 2a

figure;
set(gcf, 'color','w');
%add scatter border for evidenced solutions
%scatter(J2_HPP(IDs), J1_sed(IDs), 200,  'filled', 'MarkerEdgeColor',[ 0 0 1], 'MarkerFaceColor', [1 1 1], 'LineWidth',5);
hold on
%plot pareto front
scatter(J2_HPP, J1_sed, 50, J3_GHGs, 'filled', 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',0.01);

xlabel('HPP (GWh/yr)'); ylabel('Sediment Load (Mt/yr)'); 
%title('OBJECTIVE SPACE | 3 obj | Actual | NFE=500000 | eps=D/250');
cb = colorbar; colormap('hot'); colormap(flipud(hot)); caxis([min(J3_GHGs), max(J3_GHGs)]);
cb.Label.String = 'GHGs emissions (kgCO^2eq/yr)';

%textscatter(J2_HPP, J1_sed,  cellstr(num2str([1:length(J2_HPP)]')),'MarkerColor' ,'none','FontSize',16);

set(gca, 'FontSize',16);
hold off

%Figure 2a - detail

figure

sz_scatter_dam = 200;

[~,index] = sortrows([dams.dist2out].'); 
dams_distance = dams(index(end:-1:1)); clear index

X = [dams_distance.X]'; Y = [dams_distance.Y]';

id_dam_show = [40 42 52 60 86 123]; %1:length(dams);

name_dam = cell(size(id_dam_show));

for i=1:length(id_dam_show)
    name_dam{i} = ['(' num2str(i) ')'];
end

mapshow(MS,'Color',	[0 0.6 0.8]);
hold on
mapshow(MS(10:outlet_node),'Color',[0 0.4470 0.7410],'LineWidth',2); %show only mekong

scatter(X(id_dam_show), Y(id_dam_show), sz_scatter_dam, 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth',0.0001);

set(gcf, 'color','w');
set(gca, 'FontSize',16);

axis off
hold off

% Figure 2b

IDs = [2 86 82 280];

sz_scatter_dam = 50;

figure
for i=1:length(IDs)
    
    ID_scenario = IDs(i);

    dam_in_scenario = Portfolios(ID_scenario,:) == 1; % displays the matrix in the standard way (with 0)
        
    subplot(1,length(IDs),i);
    
    [~,index] = sortrows([dams.FromNode].'); dams_FromNode = dams(index); dam_in_scenario = dam_in_scenario(index); clear index
    X = [dams_FromNode.X]'; Y = [dams_FromNode.Y]';
    outlet_node = 128;

    range_size = [40 70*6];
    bordersize = 100;

    mapshow(MS,'Color',	[0 0.6 0.8]);
    hold on
    mapshow(MS(10:outlet_node),'Color',[0 0.4470 0.7410],'LineWidth',2); %show only mekong

    %determine dam size based on HPP
    HPPprod = [dams_FromNode.MeanAnnual]';
    size_scatter = (HPPprod - min(HPPprod)) /(max(HPPprod) - min(HPPprod))*diff(range_size) + min(range_size);

    %plot existing dams
    index_ex = [~cellfun(@(x)strcmp(x, 'P'),{dams_FromNode.Status})]==1;
    scatter(X(index_ex), Y(index_ex), size_scatter(index_ex), 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth',0.0001);

    %plot planned dams not included in portfolio
%     index_pl = and(dam_in_scenario==0, [cellfun(@(x)strcmp(x, 'P'),{dams_FromNode.Status})]==1);
%     scatter(X(index_pl), Y(index_pl), size_scatter(index_pl), 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 1 0], 'LineWidth',0.0001);

    %plot planned and included dams in portfolio
    index_pl = and(dam_in_scenario, [cellfun(@(x)strcmp(x, 'P'),{dams_FromNode.Status})]==1);
    scatter(X(index_pl), Y(index_pl), size_scatter(index_pl), 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'LineWidth',0.0001);

    %plot border for dams on Mekong
    IDreachMekong =10:outlet_node;
    index_damonMekong = sum([dams_FromNode.FromNode]' == IDreachMekong,2)' ==1;
    %index_damonMekong = and(or(index_ex, dam_in_scenario), sum([dams_FromNode.FromNode]' == IDreachMekong,2)' ==1);
    scatter(X(index_damonMekong), Y(index_damonMekong), size_scatter(index_damonMekong)/5, 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth',0.0001);
    axis off

end

set(gcf, 'color','w');
set(gca, 'FontSize',16);

hold off

%% Figure 3 

% load data

load('Borg_inputs');
load('Borg_results_multisim.mat')

% "Actual" simulation mode -> only consider "Planned" dams
outlet = 128;
index_dams_P = [];
pr_detail = [0.1 0.5 0.9];

for iii = 1:length(dams_name)
    
    dams_name(iii).dist2out = Network.Downstream.Distance{dams_name(iii).FromNode}(outlet);
    if strcmp(dams_name(iii).Status, 'P')
        index_dams_P = [index_dams_P; iii];
    end
end

dams_included = dams_name(index_dams_P);

clear iii outlet

plot_id = 2;

if plot_id == 1
    dam_plot = 1:size(prob_dam_inclusion,1);
else
    dam_plot = [48,8,7,9,55,19,18,20,21,22,16,54,50,51,13,12,14,15,60];
end

% HP_scenario_percentage_pos contains for each considered dam (row) and percentage in pr_detail (column) the position in range where the dam falls more closely 
HP_scenario_percentage_pos = zeros(length(index_dams_P),length(pr_detail));

for d= 1:length(index_dams_P)
for s=1:length(pr_detail)
    HP_scenario_percentage_pos(d,s) = (HP_scenario_percentage(d,s) - min(range_HPP))/(max(range_HPP) - min(range_HPP))*(length(range_HPP)-1)+0.5 ;
end
end

% Figure 3a

figure

colormap(flip(gray(256)))
imagesc(prob_dam_inclusion(dam_plot,:))

cb = colorbar;

cb.Label.String = '% of scenarios in which dam X is included';
cb.Label.FontSize = 20; cb.Label.FontWeight = 'bold';

hold on
sz_mrk = 100;
 scatter(HP_scenario_percentage_pos(dam_plot,1),0.8:1:length(dam_plot),sz_mrk,[0 1 0],'filled','hexagram')
 scatter(HP_scenario_percentage_pos(dam_plot,2),1:1:length(dam_plot),sz_mrk,'r','filled','hexagram')
 scatter(HP_scenario_percentage_pos(dam_plot,3),1.2:1:length(dam_plot)+0.2,sz_mrk,[0 0 1],'filled','hexagram')

xticks(0.5:length(range_HPP)+1)
xticklabels(range_HPP'/1e5)

yticks(1:length(dam_plot))

for i=1:length(dam_plot)
    name_label{i} = ['(' num2str(i) ') ' dams_included(dam_plot(i)).Name];
end
yticklabels(name_label)

ax=gca; ax.XAxis.Exponent = 5;

xlabel('Hydropower generation [10^5 GWh/yr]');

set(gcf, 'color','w');
set(gca, 'FontSize',16);

clear ax cb D 

% Figure 3b

figure

plot_id = 2;

if plot_id == 1
    dam_plot = 1:size(prob_dam_inclusion,1);
else
    dam_plot = [48,8,7,9,55,19,18,20,21,22,16,54,50,51,13,12,14,15,60];
end

data_plot = cell(1,length(pr_detail));
for i=1:length(pr_detail)
    data_plot{i} = cell2mat(cellfun(@(x)x(i,:),HP_scenario_percentage_allsim,'UniformOutput',0));
end

hold on 

B = boxplot(flip(data_plot{3}(dam_plot,:))','symbol', '.','Orientation','horizontal','BoxStyle','filled','Notch' ,'marker','Colors',[0 0 1],'Position',1:1:length(dam_plot));

xlim([min(range_HPP) max(range_HPP)])

xticks(range_HPP)
yticklabels(flip({dams_included(dam_plot).Name}))

fdamplot = flip(dam_plot);
for i=1:length(dam_plot)
    name_label{i} = ['(' num2str(length(dam_plot)-i+1) ') ' dams_included(fdamplot(i)).Name];
end
yticklabels(name_label)

 n = findobj(B,'tag','Outliers');

for j = 1:numel(n)

   n(j).MarkerEdgeColor = [1 1 1 ];
end
 
% Figure 3c

figure;
outlet = 128;

X = [dams_included.X]'; Y = [dams_included.Y]';

dam_plot = [48,8,7,9,55,19,18,20,21,22,16,54,50,51,13,12,14,15,60];

mapshow(MS,'Color',	[0 0.6 0.8]);
hold on
mapshow(MS(10:outlet),'Color',[0 0.4470 0.7410],'LineWidth',2); %show only mekong

scatter(X(dam_plot(end-5:end)), Y(dam_plot(end-5:end)),40, 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0.7 0.7], 'LineWidth',1);

scatter(X(dam_plot(1:end-6)), Y(dam_plot(1:end-6)),140, 'v', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'LineWidth',0.0001);

set(gcf, 'color','w');
set(gca, 'FontSize',16);

axis off
hold off

%% Figure 4

% load data

load('Borg_inputs');
load('Borg_results_multisim.mat')

outlet = 128;
index_dams_P = [];
index_dams_E = [];
index_dams_P_bin = zeros(size(dams))';
index_dams_E_bin = zeros(size(dams))';

for iii = 1:length(dams)
    
    dams(iii).dist2out = Network.Downstream.Distance{dams(iii).FromNode}(outlet);
    if strcmp(dams(iii).Status, 'P')
        index_dams_P = [index_dams_P; iii];
        index_dams_P_bin(iii) = 1;
    else
        index_dams_E = [index_dams_E; iii];
        index_dams_E_bin(iii) = 1;
    end
end
clear iii

dams_included = dams(index_dams_P);
index_dams_P_bin = index_dams_P_bin ==1;
index_dams_E_bin = index_dams_E_bin ==1;

ObjectID = [dams.OBJECTID]'; X = [dams.X]'; Y = [dams.Y]';

data_plot = cell(1,length(pr_detail));
for i=1:length(pr_detail)
    data_plot{i} = cell2mat(cellfun(@(x)x(i,:),HP_scenario_percentage_allsim,'UniformOutput',0));
end

clear prob_9_prc
for d = 1:length(index_dams_P)
    prob_50_prc(d,:) = [ median(data_plot{3}(d,:)), std(data_plot{3}(d,:)) ];
end

%classify dams via color
%cMap = [interp1([0;1],[0 0.3 0.6; 0.8 0.8 0 ],linspace(0,1,7)); interp1([0;1],[0.8 0.8 0; 1 0 0],linspace(0,1,7))] ;
cMap = [interp1([0;1],[ 0.8 0.8 0 ;0 0.4 0.4],linspace(0,1,7)); interp1([0;1],[0 0.4 0.4; 0.4 0 1],linspace(0,1,7))] ;
mean_value_range = range_HPP(1:end-1)+mean(diff(range_HPP))/2;

[~,ix] = min(abs(mean_value_range-prob_50_prc(:,1)),[],2);
cMap_order = cMap(ix,:); clear ix mean_value_range

% classify dams via uncertanty range
sz_scatter_dam = 70;
range_size = [sz_scatter_dam*1.2 sz_scatter_dam*4];
size_scatter = (prob_50_prc(:,2) - min(prob_50_prc(:,2))) /(max(prob_50_prc(:,2)) - min(prob_50_prc(:,2)))*diff(range_size) + min(range_size);

%plot map
figure;

mapshow(MS,'Color',	[0 0.6 0.8]);
hold on
mapshow(MS(10:outlet),'Color',[0 0.4470 0.7410],'LineWidth',2); %show only mekong

%plot dams
range_HPP_dams = [2000:4000:12000];
[~,ix] = min(abs([dams.MeanAnnual]'-range_HPP_dams),[],2);
[~,ix_pl] = min(abs([dams_included.MeanAnnual]'-range_HPP_dams),[],2);
markers = {'o','^','d'};
sz_scatter_dam = [sz_scatter_dam*0.8 sz_scatter_dam*1.2 sz_scatter_dam*1.2];

for i=1:length(range_HPP_dams) 
    
    index_pl_full = [~cellfun(@(x)strcmp(x, 'P'),{dams.Status})'.* (ix == i)]==1;
    scatter(X(index_pl_full), Y(index_pl_full), sz_scatter_dam(i)/3, markers{i}, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth',0.0001);
    
end
    
for i=1:length(range_HPP_dams) 
    
    index_pl_full = [cellfun(@(x)strcmp(x, 'P'),{dams.Status})'.* (ix == i)]==1;
    index_pl_incuded = ix_pl == i;
    scatter(X(index_pl_full), Y(index_pl_full), size_scatter(index_pl_incuded),'filled','Marker',markers{i},'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.0001);
    scatter(X(index_pl_full), Y(index_pl_full), sz_scatter_dam(i), cMap_order(index_pl_incuded,:),'filled','Marker',markers{i}, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.05);
    
end

bfr = +1e5;
xlim([min([MS.x_FN])-bfr max([MS.x_TN])+bfr])
ylim([min([MS.y_TN])-bfr max([MS.y_FN])+bfr])

%plot colorbar
colormap(cMap)
cb = colorbar;

cb.Label.String = 'Lowest HPP production in which dam X is included for 90% of the scenarios [10^5 GWh/yr]';
cb.Label.FontSize = 16; cb.Label.FontWeight = 'bold';
cb.Ticks = linspace(0, 1,length(range_HPP)) ;
cb.TickLabels = num2cell(range_HPP/1e5) ; 
 
clear scatter_1; clear scatter_2; clear scatter_3; clear bfr;

set(gcf, 'color','w');
set(gca, 'FontSize',16);

axis off
hold off
