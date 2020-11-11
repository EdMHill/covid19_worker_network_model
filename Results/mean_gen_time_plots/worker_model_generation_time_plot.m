% PURPOSE:
% Check Rt and mean generation times with no intervention measures in use
% for the worker model
%--------------------------------------------------------------------------

clear variables

%% Load data

load('../worker_model_output_RNGseed_svty_#1.mat','mean_init_generation_time_save')

% Put into 2D array. Row by replicate, column by RNGseed
mean_init_gen_time_data = squeeze(mean_init_generation_time_save);

% Set titles for the plots
title_vec = {'RNGseed: 100','RNGseed: 200','RNGseed: 300'};

%% Get mean generation time statistics

% Compute mean and percentiles
mean_of_mean_init_gen_time = mean(mean_init_gen_time_data,1);

n_configs = numel(mean_of_mean_init_gen_time);
prct_mean_init_gen_time_data = zeros(3,n_configs);
for config_itr = 1:n_configs
    prct_mean_init_gen_time_data(:,config_itr) = prctile(mean_init_gen_time_data(:,config_itr),[2.5 50 97.5]);
end
%%
% Plot histogram
position = [10, 10, 3.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])
hold on

% Specify bin edges
edges = 3.5:0.5:10.5;

% Generate the histogram
for config_itr = 1:n_configs
    % Establish the subplot
    subplot(1,n_configs,config_itr)
    hold on
    
    % Add the histogram
    histogram(mean_init_gen_time_data(:,config_itr),edges,...
                        'LineWidth',1.5)
%     histogram(mean_init_gen_time_data(:,:,config_itr),edges,...
%                         'LineWidth',1.5)

    % Set axes labels
    xlabel('Mean generation time (days)')
    if config_itr == 1
        ylabel('Frequency')
    end
    
    % Add a title
    title(title_vec{config_itr})
    
    % Add a line for the median
    med_val = median(mean_init_gen_time_data(:,config_itr));
    %med_val = median(mean_init_gen_time_data(:,:,config_itr));
    plot([med_val med_val],[0 30],'-','Color',[0.8 0. 0.],'LineWidth',1.5)
    
    % Add a text box with value of median line
    txt = num2str(med_val,'%.2f');
    text(med_val+0.5,24,txt,'Color',[0.8 0. 0.],'Fontsize',26)
    
    % Set y limit
    ylim([0 25])

    %Specify general axis properties
    set(gca,'FontSize',26)
    set(gca,'LineWidth',1)
    box on
end

% Save figure to file
export_fig('mean_generation_time_histogram','-pdf','-r1200')