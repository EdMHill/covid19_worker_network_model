% PURPOSE:
% Plot workplace distributions
% Histogram of number of workplaces (per sector)
% Couple of specific sectors chosen and workplace sizes plotted
%--------------------------------------------------------------------------

clear

% Load data
workplace_size_data = load('network_workplace_size_dist_data.mat');

% Specify global variables
selected_replicate_ID = 1;
plot_fontsize = 26;
n_simns = 100;

%% Get number of workplace data

% n_workplace_per_sector array with row per replicate, column per work
% sector
% Take mean over columns, to give mean number of workplaces per work sector (first dimension)
avg_n_workplaces_per_sector = mean(workplace_size_data.n_workplace_per_sector,1);
single_replicate_n_workplaces_per_sector = workplace_size_data.n_workplace_per_sector(selected_replicate_ID,:);

% Set histogram bins
histogram_bins_n_workplaces = [0:5:200];

% Construct averaged data histograms
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

histogram(avg_n_workplaces_per_sector,histogram_bins_n_workplaces,'Linewidth',1.5)
xlabel('Number of workplaces (per work sector)')
ylabel('Frequency')
title('Averaged over 100 simulations')

%Specify general axis properties
ylim([0 8])
set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('avg_n_workplaces_in_sector','-pdf','-r1200') % Save the figure

% Construct single replicate histograms
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

histogram(single_replicate_n_workplaces_per_sector,histogram_bins_n_workplaces,'Linewidth',1.5)
xlabel('Number of workplaces (per work sector)')
ylabel('Frequency')
title('Single simulation replicate')

%Specify general axis properties
ylim([0 8])
set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('single_simn_n_workplaces_in_sector','-pdf','-r1200') % Save the figure

%% Workplace sizes in a couple of example sectors

% Chosen sector ID to analyse
sector_ID1 = 1;
sector_ID2 = 41;

% Set histogram bins
histogram_bins_workplace_size = [0:1:50 100 150 200];

% Iterate over each simulation
% Append to vector
workplace_size_vec_1 = [];
workplace_size_vec_2 = [];
for simn_itr = 1:n_simns
    workplace_size_vec_1 = [workplace_size_vec_1; workplace_size_data.output_workplace_sizes{simn_itr}{sector_ID1}];
    workplace_size_vec_2 = [workplace_size_vec_2; workplace_size_data.output_workplace_sizes{simn_itr}{sector_ID2}];
end

% Get probability distribtuion for each simulation
% Store in array
pmf_array_1 = zeros(n_simns,numel(histogram_bins_workplace_size)-1);
pmf_array_2 = zeros(n_simns,numel(histogram_bins_workplace_size)-1);
for simn_itr = 1:n_simns
    pmf_array_1(simn_itr,:) = histcounts(workplace_size_data.output_workplace_sizes{simn_itr}{sector_ID1},...
                                histogram_bins_workplace_size,...
                                'Normalization','probability');
    pmf_array_2(simn_itr,:) = histcounts(workplace_size_data.output_workplace_sizes{simn_itr}{sector_ID2},...
                                histogram_bins_workplace_size,...
                                'Normalization','probability');                        
end

% Take mean of pmf across all simulation replicates
mean_pmf_1 = mean(pmf_array_1,1);
mean_pmf_2 = mean(pmf_array_2,1);

%Normalise
normalised_pmf_1 = mean_pmf_1./sum(mean_pmf_1);
normalised_pmf_2 = mean_pmf_2./sum(mean_pmf_2);

% Get 50+ values
bar_pmf_1 = [normalised_pmf_1(1:end-3) sum(normalised_pmf_1(end-2:end))];
bar_pmf_2 = [normalised_pmf_2(1:end-3) sum(normalised_pmf_2(end-2:end))];

%%
% Construct bar plots averaged over replicates (as probability density)
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

% Bar x-vals setup
bar_x_vals = 0:1:50;

% Set y limits
ylims_1 = [0 0.3];
ylims_2 = [0 0.45];

%histogram(normalised_pmf_1,histogram_bins_workplace_size,...
%            'Normalization','probability',...
%            'Linewidth',1.5)
bar(bar_x_vals,bar_pmf_1,...
            'Linewidth',1.5,...
            'FaceColor',[0.4000, 0.6667, 0.8431])       
xlabel('Workplace sizes (Agricultural sector)')
ylabel('Probability')
xticks([0 10 20 30 40 50])
xticklabels({'0','10','20','30','40','50+'})
ylim(ylims_1)

%Specify general axis properties
set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('avg_pmf_workplace_size_sectorID1','-pdf','-r1200') % Save the figure

position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

% histogram(normalised_pmf_,histogram_bins_workplace_size,...
%             'Normalization','probability',...
%             'Linewidth',1.5)
bar(bar_x_vals,bar_pmf_2,...
            'Linewidth',1.5,...
            'FaceColor',[0.4000, 0.6667, 0.8431]) 
xlabel('Workplace sizes (Personal Services sector)')
ylabel('Probability')
xticks([0 10 20 30 40 50])
xticklabels({'0','10','20','30','40','50+'})
ylim(ylims_2)

%Specify general axis properties
set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('avg_pmf_workplace_size_sectorID41','-pdf','-r1200') % Save the figure

%%
% Construct bar plot for single replicate (as probability density)
% For sector_ID1
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

% histogram(workplace_size_data.output_workplace_sizes{selected_replicate_ID}{sector_ID1},...
%              histogram_bins_workplace_size,...
%             'Normalization','probability',...
%             'Linewidth',1.5)
counts_1 = histcounts(workplace_size_data.output_workplace_sizes{selected_replicate_ID}{sector_ID1},[bar_x_vals 200]);
pmf_single_rep_1 = counts_1./sum(counts_1);
bar(bar_x_vals,pmf_single_rep_1,...
            'Linewidth',1.5,...
            'FaceColor',[0.4000, 0.6667, 0.8431])    
xlabel('Workplace sizes (Agricultural sector)')
ylabel('Probability')
xticks([0 10 20 30 40 50])
xticklabels({'0','10','20','30','40','50+'})
ylim(ylims_1)

set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('one_simn_pmf_workplace_size_sectorID1','-pdf','-r1200') % Save the figure

% Construct histogram for single replicate (as probability density)
% For sector_ID2
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])

% histogram(workplace_size_data.output_workplace_sizes{selected_replicate_ID}{sector_ID2},...
%             histogram_bins_workplace_size,...
%             'Normalization','probability',...
%             'Linewidth',1.5)
counts_2 = histcounts(workplace_size_data.output_workplace_sizes{selected_replicate_ID}{sector_ID2},[bar_x_vals 200]);
pmf_single_rep_2 = counts_2./sum(counts_2);
bar(bar_x_vals,pmf_single_rep_2,...
            'Linewidth',1.5,...
            'FaceColor',[0.4000, 0.6667, 0.8431]) 
xlabel('Workplace sizes (Personal Services sector)')
ylabel('Probability')
xticks([0 10 20 30 40 50])
xticklabels({'0','10','20','30','40','50+'})
ylim(ylims_2)

set(gca,'FontSize',plot_fontsize)
set(gca,'LineWidth',1)
box on

export_fig('one_simn_pmf_workplace_size_sectorID41','-pdf','-r1200') % Save the figure
