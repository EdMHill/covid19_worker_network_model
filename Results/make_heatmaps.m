% PURPOSE:
% Produce heatmaps to show variability in outputs across range of values for
% two input variables.
% Used with covid-secure scenario outputs
%--------------------------------------------------------------------------

function make_heatmaps(dataset,variable_names)


%% Load data
if strcmp(dataset,'COVID_secure')==1
    CS_data = load('worker_model_output_CS_intervention_combined.mat');
elseif strcmp(dataset,'COVID_secure_no_isol')==1
    CS_data = load('worker_model_output_CS_intervention_no_isol_combined.mat');
end

%% Specify global params

% Specify network size
cmax = 10000;

% Number of simulations replicates performed
n_simns = 1000;

%% Specify global plot variable values

% Set axis labels
xaxis_label = 'Maximum work team size';
yaxis_label = 'Work contact transmission risk scaling';

% Set plot fontsize
plot_fontsize = 30;

%% Specify save file identifier
if strcmp(dataset,'COVID_secure')==1
    save_file_addon = '_with_isol';
elseif strcmp(dataset,'COVID_secure_no_isol')==1
    save_file_addon = '_no_isol';
end

%% Compute desired summary statistics

if sum(strcmp(variable_names,'final_size') == 1)

    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(CS_data.numinf_combined(end,:,:));

    % Propn of infections over duration of outbreak
    final_size_propn = final_size_absolute/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(final_size_propn>0.5)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
    
    % Call plot function
    panel_title = 'Cumulative infectious case proportion > 0.5';
    add_colourbar = true;
    save_filename = ['covid_secure_plots/final_size_heatmap',save_file_addon];
    plot_heatmap(heatmap_array,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
end

if sum(strcmp(variable_names,'peak_inf') == 1)

    % Number of infections over duration of outbreak
    peak_inf_CS_data = squeeze(max(CS_data.newinf_combined))/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(peak_inf_CS_data>0.01)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
    
    % Call plot function
    panel_title = 'Peak in infectious prevalence > 0.01';
    add_colourbar = true;
    save_filename = ['covid_secure_plots/peak_inf_heatmap',save_file_addon];
    plot_heatmap(heatmap_array,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
end

if sum(strcmp(variable_names,'avg_isolation') == 1)
        
    % Time in isolation (average per person)
    offset_idx = 31;
    total_isol_CS_data = squeeze(sum(CS_data.num_isolating_combined(offset_idx:end,:,:)))/(cmax*365);

    % Check against the threshold criteria and store
    heatmap_vals = sum(total_isol_CS_data>0.01)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
    
    % Call plot function
    panel_title = 'Proportion of time in isolation > 0.01';
    add_colourbar = true;
    save_filename = ['covid_secure_plots/avg_isol_heatmap',save_file_addon];
    plot_heatmap(heatmap_array,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
end

if sum(strcmp(variable_names,'peak_isolation') == 1)

    % Peak isolations
    peak_isol_CS_data = squeeze(max(CS_data.num_isolating_combined))/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(peak_isol_CS_data>0.05)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);

    % Call plot function
    panel_title = 'Peak in proportion isolated > 0.05';
    add_colourbar = true;
    save_filename = ['covid_secure_plots/peak_isol_heatmap',save_file_addon];
    plot_heatmap(heatmap_array,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
end

if sum(strcmp(variable_names,'duration') == 1)

    % Peak isolations
    prev_data = CS_data.prevpresymp_combined + CS_data.prevasymp_combined + CS_data.prevsymp_combined;
    duration_data = zeros(size(squeeze(prev_data(1,:,:))));
    for i = 1:length(CS_data.numinf_combined(1,:,1))
        for j = 1:length(CS_data.numinf_combined(1,1,:))
            duration_data(i,j) = find(prev_data(:,i,j)>0,1,'last') - 1;        
        end
    end
    
    % Check against the threshold criteria and store
    heatmap_vals = sum(duration_data>150)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);

    % Call plot function
    panel_title = 'Outbreak duration > 150 days';
    add_colourbar = true;
    save_filename = ['covid_secure_plots/durations_heatmap',save_file_addon];
    plot_heatmap(heatmap_array,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
end
end

%% Create the heatmaps

function plot_heatmap(heatmap_data,yaxis_label,xaxis_label,panel_title,...
                            add_colourbar,plot_fontsize,save_filename)
    % Set up plot
    position = [10, 10, 2.0*550, 2.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])

    % Produce heatmap
    imagesc(heatmap_data)

    % Set y-axis labels
    ylabel(yaxis_label)

    % Set x-axis labels
    xlabel(xaxis_label)

    % % Set axis limits
    % xlim(xlim_vals)
    % ylim(ylim_vals)

    % Set axis ticks
    yticks(1:1:4)
    yticklabels({'0.25','0.50','0.75','1.00'})
    xticks(1:1:3)
    xticklabels({'2','5','10'})

    % Add title
    title(panel_title)

    % Add colourbar (if applicable)
    caxis([0 1])
    if add_colourbar == true
        c = colorbar;
        c.Label.String = 'Fraction of simulations satisfying the criteria';
        c.LineWidth = 1.5;
        c.FontSize = floor(plot_fontsize);
        set(c,'YTickLabel',{'0.0';'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9';'1.0'})
    end

    % Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on
    
    % Save the figure
    export_fig(save_filename,'-pdf','-transparent','-r1200')
end

