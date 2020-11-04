% PURPOSE:
% Plot profiles of fraction of simulations that satisfy given criteria
%--------------------------------------------------------------------------

function make_threshold_event_plots(dataset,...
                                    variable_names,...
                                    criteria_thresholds,...
                                    display_legend)
% Inputs: 
% dataset - Data from the batch of scenarios being considered (e.g. 'adherence')
% variable_names - String vector specifying what measures are being looked at
% criteria_thresholds - The condition to be realised for that outcome to be logged as exceeding the given threshold  
% display_legend/ - Whether to display the legend

% Input variable check. The number of variables being analysed and threshold criteria
% must match
if numel(variable_names) ~= numel(criteria_thresholds)
    error('The number of entries in variable_names (%d) does not match the number of entires in criteria_thresholds (%d)',...
                numel(variable_names),numel(criteria_thresholds));
end

if numel(variable_names) > 5
    error('More than five variables names have been provided. Maximum number allowed is five.')
end

%% Load data
if strcmp(dataset,'adherence')==1
    % Adherence
    input_data = load('worker_model_output_adherence_intervention_combined.mat');
    
    % Set if will be used for collective plotting with other datasets
    multi_panel_plot_flag = false;
elseif strcmp(dataset,'workpercent')==1
    % Work percent
    input_data = load('worker_model_output_workpercent_intervention_combined.mat');
    
    % Set if will be used for collective plotting with other datasets
    multi_panel_plot_flag = false;
elseif strcmp(dataset,'backwardsCT')==1
    % Backwards CT
    input_data = load('worker_model_output_amount_backwards_CT_intervention_combined.mat');
    
    % Set if will be used for collective plotting with other datasets
    multi_panel_plot_flag = false;
elseif strcmp(dataset,'workpercent_and_backwardsCT')==1
    %Synchronised and asynchronised worker patterns
    input_data = load('worker_model_output_workpercent_intervention_combined.mat');
    input_data_2 = load('worker_model_output_amount_backwards_CT_intervention_combined.mat');
    
    % Set if will be used for collective plotting with other datasets
    multi_panel_plot_flag = true;
    
    % Set up panel titles
    panel_title = {'Working from home','Infector identified'};
elseif strcmp(dataset,'worker_patterns')==1
    %Synchronised and asynchronised worker patterns
    input_data = load('worker_model_output_synchronised_changedays_intervention_combined.mat');
    input_data_2 = load('worker_model_output_variable_changedays_intervention_combined.mat');
    
    % Set if will be used for collective plotting with other datasets
    multi_panel_plot_flag = true;
    
    % Set up panel titles
    panel_title = {'Synchronised','Asynchronised'};
end

%% Specify global params

% Specify network size
cmax = 10000;

% Number of simulations replicates performed
n_simns = 1000;

% Set up trace colours for plots
colour_array = [0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.4940    0.1840    0.5560;
                0.5 0.5 0.5;
                0.4660    0.6740    0.1880;
                0.9290    0.6940    0.1250];

% Set plot fontsize
plot_fontsize = 26;

%% Initialise threshold event variables

% number of scenarios under consideration (e.g. adherence values tested)
n_scens = size(input_data.numinf_combined,3);

% Number of variables being considered
n_vars = numel(variable_names);

% Initialise the storage array
% Intention to plot each column (for criteria being assessed) 
% against the scenario configs.
threshold_vals = zeros(n_scens,n_vars);

% If plotting multiple panels collectively, set up additional storage
% arrays
if multi_panel_plot_flag == true
    n_scens_2 = size(input_data_2.numinf_combined,3);
    threshold_vals_2 = zeros(n_scens_2,n_vars);
end

% Initialise legend label string vector
legend_labels = cell(n_vars,1);

%% Compute desired summary statistics

% Find if variable is in list
if sum(strcmp(variable_names,'final_size') == 1)
    
    % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'final_size') == 1);
    
    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(input_data.numinf_combined(end,:,:));
    
    % Propn of infections over duration of outbreak
    final_size_propn = final_size_absolute/cmax;

    % Check against the threshold criteria and store
    threshold_vals(:,threshold_idx) = sum(final_size_propn>criteria_thresholds(threshold_idx))/n_simns;

    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = ['Propn infectious > ',num2str(criteria_thresholds(threshold_idx))];
    
    % If working with multiple sets of data, get similar results for those
    % sets of data
    if multi_panel_plot_flag == true
        final_size_absolute_2 = squeeze(input_data_2.numinf_combined(end,:,:));
        final_size_propn_2 = final_size_absolute_2/cmax;
        threshold_vals_2(:,threshold_idx) = sum(final_size_propn_2>criteria_thresholds(threshold_idx))/n_simns;
    end
end

if sum(strcmp(variable_names,'peak_inf') == 1)
    
    % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'peak_inf') == 1);
        
    % Peak proportion new infections
    peak_inf_input_data = squeeze(max(input_data.newinf_combined))/cmax;

    % Check against the threshold criteria and store
    threshold_vals(:,threshold_idx) = sum(peak_inf_input_data>criteria_thresholds(threshold_idx))/n_simns;

    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = ['Peak infectious > ',num2str(criteria_thresholds(threshold_idx))];
    
    % If working with multiple sets of data, get similar results for those
    % sets of data
    if multi_panel_plot_flag == true
        peak_inf_input_data_2 = squeeze(max(input_data_2.newinf_combined))/cmax;
        threshold_vals_2(:,threshold_idx) = sum(peak_inf_input_data_2>criteria_thresholds(threshold_idx))/n_simns;
    end
end

if sum(strcmp(variable_names,'inf_size_or_peak') == 1)
    
    % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'inf_size_or_peak') == 1);
        
    % Propn of infections over duration of outbreak
    final_size_absolute = squeeze(input_data.numinf_combined(end,:,:));
    final_size_propn = final_size_absolute/cmax; 
    
    % Peak proportion new infections
    peak_inf_input_data = squeeze(max(input_data.newinf_combined))/cmax;

    % Check against the threshold criteria and store
    propn_inf_indicator_array = final_size_propn>0.2;
    peak_indicator_array = peak_inf_input_data>0.01;
    conditions_satisfied = propn_inf_indicator_array + peak_indicator_array;
    threshold_vals(:,threshold_idx) = sum(conditions_satisfied>0)/n_simns;
 
    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = 'Propn infectious > 0.2 OR Peak infectious > 0.01';
end

if sum(strcmp(variable_names,'avg_isolation') == 1)
     % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'avg_isolation') == 1);
        
    % Time in isolation (average per person)
    offset_idx = 31;
    total_isol_input_data = squeeze(sum(input_data.num_isolating_combined(offset_idx:end,:,:)))/(cmax*365);

    % Check against the threshold criteria and store
    threshold_vals(:,threshold_idx) = sum(total_isol_input_data>criteria_thresholds(threshold_idx))/n_simns;

    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = ['Avg time in isolation > ',num2str(criteria_thresholds(threshold_idx))];
    
    % If working with multiple sets of data, get similar results for those
    % sets of data
    if multi_panel_plot_flag == true
        total_isol_input_data_2 = squeeze(sum(input_data_2.num_isolating_combined(offset_idx:end,:,:)))/(cmax*365);
        threshold_vals_2(:,threshold_idx) = sum(total_isol_input_data_2>criteria_thresholds(threshold_idx))/n_simns;
    end
end

if sum(strcmp(variable_names,'peak_isolation') == 1)
     % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'peak_isolation') == 1);
        
    % Peak isolations
    peak_isol_input_data = squeeze(max(input_data.num_isolating_combined))/cmax;

    % Check against the threshold criteria and store
    threshold_vals(:,threshold_idx) = sum(peak_isol_input_data>criteria_thresholds(threshold_idx))/n_simns;

    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = ['Peak propn isolated > ',num2str(criteria_thresholds(threshold_idx))];
    
    % If working with multiple sets of data, get similar results for those
    % sets of data
    if multi_panel_plot_flag == true
        peak_isol_input_data_2 = squeeze(max(input_data_2.num_isolating_combined))/cmax;
        threshold_vals_2(:,threshold_idx) = sum(peak_isol_input_data_2>criteria_thresholds(threshold_idx))/n_simns;
    end
end

if sum(strcmp(variable_names,'duration') == 1)
     % Get the placement index of the variable name from the vector variable_names
    threshold_idx = find(strcmp(variable_names,'duration') == 1);
    
    % Peak isolations
    prev_data = input_data.prevpresymp_combined + input_data.prevasymp_combined + input_data.prevsymp_combined;
    duration_data = zeros(size(squeeze(prev_data(1,:,:))));
    for i = 1:length(input_data.numinf_combined(1,:,1))
        for j = 1:length(input_data.numinf_combined(1,1,:))
            duration_data(i,j) = find(prev_data(:,i,j)>0,1,'last');        
        end
    end
    
    % Check against the threshold criteria and store
    threshold_vals(:,threshold_idx) = sum(duration_data>criteria_thresholds(threshold_idx))/n_simns;

    % Assign legend label to plot associated vector
    legend_labels{threshold_idx} = ['Outbreak duration > ',num2str(criteria_thresholds(threshold_idx)),' days'];
    
    % If working with multiple sets of data, get similar results for those
    % sets of data
    if multi_panel_plot_flag == true
        prev_data_2 = input_data_2.prevpresymp_combined + input_data_2.prevasymp_combined + input_data_2.prevsymp_combined;
        duration_data_2 = zeros(size(squeeze(prev_data_2(1,:,:))));
        for i = 1:length(input_data_2.numinf_combined(1,:,1))
            for j = 1:length(input_data_2.numinf_combined(1,1,:))
                duration_data_2(i,j) = find(prev_data_2(:,i,j)>0,1,'last');        
            end
        end
        threshold_vals_2(:,threshold_idx) = sum(duration_data_2>criteria_thresholds(threshold_idx))/n_simns;
    end
end

%% Adherence sensitivity
if strcmp(dataset,'adherence')==1
       
    % Set x-axis label
    xaxis_label = 'Adherence probability';
    
    % Set up xticks
    xticks_vals = 0:1:10;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
    
    % Set up xaxis limits
    xlim_vals = [-0.2 10.2];
    ylim_vals = [0 1.02];
    
    % Set filename prefix
    save_filename = 'threshold_event_plots/threshold_events_adherence_scen';
elseif strcmp(dataset,'workpercent')==1
    % Set x-axis label
    xaxis_label = 'Work from home probability';
    
    % Set up xticks
    xticks_vals = 0:1:11;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','N-U'};
    
    % Set up xaxis limits
    xlim_vals = [-0.2 11.2];
    ylim_vals = [0 1.02];
    
    % Set filename prefix
    save_filename = 'threshold_event_plots/threshold_events_workpercent_scen';
elseif strcmp(dataset,'backwardsCT')==1
    % Set x-axis label
    xaxis_label = 'Backward contact tracing: Infector identified probability';
    
    % Set up xticks
    xticks_vals = 0:1:10;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
    
    % Set up xaxis limits
    xlim_vals = [-0.2 10.2];
    ylim_vals = [0 1.02];
    
    % Set filename prefix
    save_filename = 'threshold_event_plots/threshold_events_backwardsCT_scen';
elseif strcmp(dataset,'workpercent_and_backwardsCT')==1
    % Set x-axis label
    xaxis_label = 'Work from home probability';
    
    % Set up xticks
    xticks_vals = 0:1:11;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','N-U'};
    xticks_vals_2 = 0:1:10;
    xticks_labels_2 = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
    
    
    % Set up xaxis limits
    xlim_vals = [-0.2 11.2];
    xlim_vals_2 = [-0.2 10.2];
    ylim_vals = [0 1.02];
    
    % Set filename prefix
    save_filename = 'threshold_event_plots/threshold_events_workpercent_and_backwardsCT_scen';
elseif strcmp(dataset,'worker_patterns')==1
    % Set x-axis label
    xaxis_label = 'Days per week at the workplace';
    
    % Set up xticks
    xticks_vals = 0:1:5;
    xticks_labels = {'0','1','2','3','4','5'};
    xticks_vals_2 = 0:1:5;
    xticks_labels_2 = {'0','1','2','3','4','5'};
    
    % Set up xaxis limits
    xlim_vals = [-0.2 5.2];
    xlim_vals_2 = [-0.2 5.2];
    ylim_vals = [0 max(threshold_vals(:))*1.05];
    
    % Set filename prefix
    save_filename = 'threshold_event_plots/threshold_events_worker_pattern_scen';
end

%% Generate the plot

if ispc==1
    scale = 1/0.75;
else
    scale = 1;
end

% For some plots, have double panel, with space for legend to the left
if multi_panel_plot_flag == true
    % Set up plot
    position = [10, 10, 3.5*550*scale, 1.5*450*scale];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    
    % Plot synchronised worker pattern data
    subplot(1,10,1:4)
    hold on
    
    % Plot the traces
    for trace_itr = 1:n_vars
        if trace_itr == 1
            plot(xticks_vals,threshold_vals(:,trace_itr),'x-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 2
            plot(xticks_vals,threshold_vals(:,trace_itr),'v--',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 3
            plot(xticks_vals,threshold_vals(:,trace_itr),'o:',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
         elseif trace_itr == 4
            plot(xticks_vals,threshold_vals(:,trace_itr),'*-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 5
            plot(xticks_vals,threshold_vals(:,trace_itr),'s-.',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        end
    end

    % Set y-axis labels
    ylabel('Proportion of simulations')

    % Set x-axis labels
    xlabel(xaxis_label)

    % Set axis limits
    xlim(xlim_vals)
    ylim(ylim_vals)

    % Set axis ticks
    yticks(0:0.1:1)
    xticks(xticks_vals)
    xticklabels(xticks_labels)
    
    % Set panel title
    title(panel_title{1})

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on
    
    % Plot asynchronised worker pattern data
    subplot(1,10,5:8)
    hold on
    
    % Plot the traces
    for trace_itr = 1:n_vars
        if trace_itr == 1
            plot(xticks_vals_2,threshold_vals_2(:,trace_itr),'x-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 2
            plot(xticks_vals_2,threshold_vals_2(:,trace_itr),'v--',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 3
            plot(xticks_vals_2,threshold_vals_2(:,trace_itr),'o:',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 4
            plot(xticks_vals_2,threshold_vals_2(:,trace_itr),'*-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 5
            plot(xticks_vals_2,threshold_vals_2(:,trace_itr),'s-.',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        end
    end

    % Set axis labels
    xlabel(xaxis_label)
    set(gca,'YTickLabel',[]);

    % Set axis limits
    xlim(xlim_vals_2)
    ylim(ylim_vals)

    % Set axis ticks
    xticks(xticks_vals_2)
    xticklabels(xticks_labels_2)
    
    % Set panel title
    title(panel_title{2})

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on
    
    % Add legend, if applicable
    if display_legend == true
        H = gobjects(n_vars,1);
        for data_itr = 1:n_vars
            if data_itr == 1
                H(data_itr) = plot(NaN,NaN,'x-','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 2
                H(data_itr) = plot(NaN,NaN,'v--','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 3
               H(data_itr) = plot(NaN,NaN,'o:','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 4
               H(data_itr) = plot(NaN,NaN,'*-','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 5
               H(data_itr) = plot(NaN,NaN,'s-.','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            end
        end
        legend(H,...
            'LineWidth',1.5,...
            'FontSize',floor(plot_fontsize)-2,...
            'Position',[0.76 0.43 0.172 0.177]);
    end
    
else
    % Set up plot
    position = [10, 10, 2.5*550*scale, 1.5*450*scale];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    % Plot the traces
    for trace_itr = 1:n_vars
        if trace_itr == 1
            plot(xticks_vals,threshold_vals(:,trace_itr),'x-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 2
            plot(xticks_vals,threshold_vals(:,trace_itr),'v--',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 3
            plot(xticks_vals,threshold_vals(:,trace_itr),'o:',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
       elseif trace_itr == 4
            plot(xticks_vals,threshold_vals(:,trace_itr),'*-',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        elseif trace_itr == 5
            plot(xticks_vals,threshold_vals(:,trace_itr),'s-.',...
                'Color',colour_array(trace_itr,:),...
                'MarkerSize',12,...
                'LineWidth',1.5)
        end
    end

    % Set y-axis labels
    ylabel('Proportion of simulations')

    % Set x-axis labels
    xlabel(xaxis_label)

    % Set axis limits
    xlim(xlim_vals)
    ylim(ylim_vals)

    % Set axis ticks
    yticks(0:0.1:1)
    xticks(xticks_vals)
    xticklabels(xticks_labels)

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on

    % Add legend, if applicable
    if display_legend == true
        H = gobjects(n_vars,1);
        for data_itr = 1:n_vars
            if data_itr == 1
                H(data_itr) = plot(NaN,NaN,'x-','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 2
                H(data_itr) = plot(NaN,NaN,'v--','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 3
               H(data_itr) = plot(NaN,NaN,'o:','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 4
               H(data_itr) = plot(NaN,NaN,'*-','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            elseif data_itr == 5
               H(data_itr) = plot(NaN,NaN,'s-.','Color',colour_array(data_itr,:),...
                'DisplayName',legend_labels{data_itr},...
                'MarkerSize',12,...
                'LineWidth',1.5);
            end
        end
        legend(H,...
            'LineWidth',1.5,...
            'FontSize',floor(plot_fontsize)-2,...
            'Location','eastoutside');
    end
end

% Save figure to file
export_fig(save_filename,'-pdf','-transparent','-r1200')
end