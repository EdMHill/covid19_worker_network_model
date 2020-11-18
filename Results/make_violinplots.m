% PURPOSE:
% Make violin plots
%--------------------------------------------------------------------------

% clear variables

function make_violinplots(dataset,variablename)

%% ADD PATH DEPENDENCIES
addpath('../../Matlab/Violinplot-Matlab')

%% Load data
% Adherence
adherence_data = load('worker_model_output_adherence_intervention_combined.mat');

% Work percent
workpercent_data = load('worker_model_output_workpercent_intervention_combined.mat');

% Backwards CT
backwards_CT_data = load('worker_model_output_amount_backwards_CT_intervention_combined.mat');

% Synchronised
synch_data = load('worker_model_output_synchronised_changedays_intervention_combined.mat');

% Non-synchronised
asynch_data = load('worker_model_output_variable_changedays_intervention_combined.mat');

%% Specify global params

% Specify network size
cmax = 10000;

%% Compute desired summary statistics
if strcmp(variablename,'final_size')
    % Number of infections from day 15 onwards (i.e. timestep 16)
    offset_idx = 16;
    adherence_final_size = squeeze(adherence_data.numinf_combined(end,:,:) - adherence_data.numinf_combined(offset_idx,:,:));
    workpercent_final_size = squeeze(workpercent_data.numinf_combined(end,:,:) - workpercent_data.numinf_combined(offset_idx,:,:));
    backwards_CT_final_size = squeeze(backwards_CT_data.numinf_combined(end,:,:) - backwards_CT_data.numinf_combined(offset_idx,:,:));
    synch_final_size = squeeze(synch_data.numinf_combined(end,:,:) - synch_data.numinf_combined(offset_idx,:,:));
    asynch_final_size = squeeze(asynch_data.numinf_combined(end,:,:) - asynch_data.numinf_combined(offset_idx,:,:));
    
    % Proportion of infections from day 15 onwards (i.e. timestep 16)
    adherence_input_data = adherence_final_size/cmax;
    workpercent_input_data = workpercent_final_size/cmax;
    backwards_CT_input_data = backwards_CT_final_size/cmax;
    synch_input_data = synch_final_size/cmax;
    asynch_input_data = asynch_final_size/cmax;
elseif strcmp(variablename,'peak_inf')
    % Peak proportion new infections
    adherence_input_data = squeeze(max(adherence_data.newinf_combined))/cmax;
    workpercent_input_data = squeeze(max(workpercent_data.newinf_combined))/cmax;
    backwards_CT_input_data = squeeze(max(backwards_CT_data.newinf_combined))/cmax;
    synch_input_data = squeeze(max(synch_data.newinf_combined))/cmax;
    asynch_input_data = squeeze(max(asynch_data.newinf_combined))/cmax;
elseif strcmp(variablename,'total_isolation')
    % total number of isolation-days from day 15 onwards (i.e. timestep 16)
    offset_idx = 16;
    adherence_input_data = squeeze(sum(adherence_data.num_isolating_combined(offset_idx:end,:,:)));
    workpercent_input_data = squeeze(sum(workpercent_data.num_isolating_combined(offset_idx:end,:,:)));
    backwards_CT_input_data = squeeze(sum(backwards_CT_data.num_isolating_combined(offset_idx:end,:,:)));
    synch_input_data = squeeze(sum(synch_data.num_isolating_combined(offset_idx:end,:,:)));
    asynch_input_data = squeeze(sum(asynch_data.num_isolating_combined(offset_idx:end,:,:)));
elseif strcmp(variablename,'peak_isolation')
    % peak number isolating
    adherence_input_data = squeeze(max(adherence_data.num_isolating_combined))/cmax;
    workpercent_input_data = squeeze(max(workpercent_data.num_isolating_combined))/cmax;
    backwards_CT_input_data = squeeze(max(backwards_CT_data.num_isolating_combined))/cmax;
    synch_input_data = squeeze(max(synch_data.num_isolating_combined))/cmax;
    asynch_input_data = squeeze(max(asynch_data.num_isolating_combined))/cmax;
elseif strcmp(variablename,'duration')
    % outbreak duration
    adherence_data.prev = adherence_data.prevpresymp_combined + adherence_data.prevasymp_combined + adherence_data.prevsymp_combined;
    workpercent_data.prev = workpercent_data.prevpresymp_combined + workpercent_data.prevasymp_combined + workpercent_data.prevsymp_combined;
    backwards_CT_data.prev = backwards_CT_data.prevpresymp_combined + backwards_CT_data.prevasymp_combined + backwards_CT_data.prevsymp_combined;
    synch_data.prev = synch_data.prevpresymp_combined + synch_data.prevasymp_combined + synch_data.prevsymp_combined;
    asynch_data.prev = asynch_data.prevpresymp_combined + asynch_data.prevasymp_combined + asynch_data.prevsymp_combined;
    adherence_input_data = zeros(size(squeeze(adherence_data.prev(1,:,:))));
    workpercent_input_data = zeros(size(squeeze(workpercent_data.prev(1,:,:))));
    backwards_CT_input_data = zeros(size(squeeze(backwards_CT_data.prev(1,:,:))));
    synch_input_data = zeros(size(squeeze(synch_data.prev(1,:,:))));
    asynch_input_data = zeros(size(squeeze(asynch_data.prev(1,:,:))));
    for i = 1:length(adherence_data.numinf_combined(1,:,1))
        for j = 1:length(adherence_data.numinf_combined(1,1,:))
%             adherence_input_data(i,j) = find(adherence_data.numinf_combined(:,i,j)>0,1,'last');
%             workpercent_input_data(i,j) = find(workpercent_data.numinf_combined(:,i,j)>0,1,'last');
%             backwards_CT_input_data(i,j) = find(backwards_CT_data.numinf_combined(:,i,j)>0,1,'last');
            adherence_input_data(i,j) = find(adherence_data.prev(:,i,j)>0,1,'last');
            backwards_CT_input_data(i,j) = find(backwards_CT_data.prev(:,i,j)>0,1,'last');
        end
    end
    for i = 1:length(workpercent_data.numinf_combined(1,:,1))
        for j = 1:length(workpercent_data.numinf_combined(1,1,:))
            workpercent_input_data(i,j) = find(workpercent_data.prev(:,i,j)>0,1,'last');
        end
    end
    for i = 1:length(synch_data.numinf_combined(1,:,1))
        for j = 1:length(synch_data.numinf_combined(1,1,:))
            synch_input_data(i,j) = find(synch_data.prev(:,i,j)>0,1,'last');
            asynch_input_data(i,j) = find(asynch_data.prev(:,i,j)>0,1,'last');
        end
    end
end

%% Adherence sensitivity
if strcmp(dataset,'adherence')==1
    % Set the scenario name
    scen_name = "sweep_adherence";
    
    % Specify input data
    input_data = {adherence_input_data};
    
    % Set plot style
    plot_style = 'violin';
    %plot_style = 'boxplot';
    
    % Set colour for violin plots
    colour_vec = [0.8500    0.3250    0.0980];
    
    % Set x-axis label
    xaxis_label = 'Adherence probability';
    
    % Set up xticks
    xticks_vals = 0:1:10;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
    
    % Set up xaxis limits
    xlim_vals = [-0.5 10.5];
    ylim_vals = [0 max(input_data{1}(:))*1.05];
    
    % Set filename prefix
    save_filename_initial = ['adherence_sens_plots/worker_model_sens_adherence_',variablename];
    
    % Set plot fontsize
    plot_fontsize = 22;
    
    % Call function to produce plots
    propn_plot_type = false;
    
    % Do not display a legend
    display_legend = false;
    legend_label = {};
elseif strcmp(dataset,'workpercent')==1
    % Set the scenario name
    scen_name = "workpercent";
    
    % Specify input data
    input_data = {workpercent_input_data};
    
    % Set plot style
    plot_style = 'violin';
    
    % Set colour for violin plots
    colour_vec = [1 0.7 0.7];
    
    % Set x-axis label
    xaxis_label = 'Proportion';
    
    % Set up xticks
    xticks_vals = 0:1:11;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','N-U'};
    
    % Set up xaxis limits
    xlim_vals = [-0.5 11.5];
    ylim_vals = [0 max(input_data{1}(:))*1.05];
    
    % Set filename prefix
    save_filename_initial = ['workpercent_plots/worker_model_sens_workpercent_',variablename];
    
    % Set plot fontsize
    plot_fontsize = 26;
    
    % Call function to produce plots
    propn_plot_type = false;
    
    % Do not display a legend
    display_legend = false;
    legend_label = {};
elseif strcmp(dataset,'workpercent_and_backwardsCT')==1
    %% Workpercent and backwards CT comparison
    
    % Set the scenario name
    scen_name = "workpercent_backwardsCT_compare";
    
    % Specify input data
    input_data = {workpercent_input_data,backwards_CT_input_data};
    
    % Set plot style
    plot_style = 'violin';
    %plot_style = 'boxplot';
    
    % Set colour for violin plots
    colour_vec = [1    0.7    0.7;
        0.3 0.3 1
        ];
    
    % Set x-axis label
    xaxis_label = 'Proportion';
    
    % Set up xticks
    xticks_vals = 0:1:11;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','N-U'};
    
    % Set up xaxis limits
    xlim_vals = [-0.5 11.5];
    ylim_vals = [0 max(input_data{1}(:))*1.05];
    
    % Set filename prefix
    save_filename_initial = ['misc_plots/worker_model_sens_workpercent_backwardsCT_',variablename];
    
    % Set plot fontsize
    plot_fontsize = 26;
    
    % Call function to produce plots
    propn_plot_type = false;
    
    % Set legend variables
    if strcmp(variablename,'final_size')
        display_legend = true;
    else
        display_legend = false;
    end
    legend_label = {'Working from home','Backward contact tracing: Infector identified'};
elseif strcmp(dataset,'worker_patterns')==1
    %% Synch & asynch worker pattern comparison
    
    % Specify input data
    input_data = {synch_input_data,asynch_input_data};
    
    % Set plot style
    plot_style = 'violin';
    %plot_style = 'boxplot';
    
    % Set colour for violin plots
    colour_vec = [0.4 0 0;
        0.2 0.8 0.8
        ];
    
    % Set x-axis label
    xaxis_label = 'Days per week at the workplace';
    
    % Set up xticks
    xticks_vals = 0:1:5;
    xticks_labels = {'0','1','2','3','4','5'};
    
    % Set up xaxis limits
    xlim_vals = [-0.5 5.5];
    ylim_vals = [0 max(input_data{1}(:))*1.05];
    
    % Set filename prefix
    save_filename_initial = ['worker_pattern_plots/worker_model_sens_worker_patterns_',variablename];
    
    % Set plot fontsize
    plot_fontsize = 26;
    
    % Call function to produce plots
    propn_plot_type = false;
    
    % Set legend variables
    if strcmp(variablename,'final_size')
        display_legend = true;
    else
        display_legend = false;
    end
    legend_label = {'Synchronous','Asynchronous'};
end

% Set up y-axis labels per statistic.
% If required, normalise the data
if propn_plot_type == true
    if strcmp(variablename,'final_size')
        yaxis_label = 'Relative proportion infected';
    elseif strcmp(variablename,'peak_inf')
        yaxis_label = 'Relative peak proportion new infections';
    elseif strcmp(variablename,'total_isolation')    
        yaxis_label = 'Relative total isolation-days';    
    elseif strcmp(variablename,'peak_isolation')
        yaxis_label = 'Relative peak proportion in isolation';
    elseif strcmp(variablename,'duration')
        yaxis_label = 'Relative outbreak duration';
    end
else
    if strcmp(variablename,'final_size')
        yaxis_label = 'Additional fraction of the network infected';
    elseif strcmp(variablename,'peak_inf')
        yaxis_label = 'Peak proportion new infections';
    elseif strcmp(variablename,'total_isolation')
        yaxis_label = 'Total isolation-days';  
    elseif strcmp(variablename,'peak_isolation')
        yaxis_label = 'Peak proportion in isolation';
    elseif strcmp(variablename,'duration')
        yaxis_label = 'Outbreak duration';
    end
end

%% Call plot function

% If required, normalise values relative to first column of data
if propn_plot_type == true
    baseline_vals = input_data(:,1);
    input_data = input_data./baseline_vals;
end

generate_sensitivity_plot(input_data,...
        propn_plot_type,...
        plot_style,...
        colour_vec,...
        yaxis_label,...
        xaxis_label,...
        xticks_vals,...
        xticks_labels,...
        ylim_vals,...
        xlim_vals,...
        display_legend,...
        legend_label,...
        save_filename_initial,...
        plot_fontsize,...
        dataset)
end

%% Function to generate the violin/boxplots plots
function generate_sensitivity_plot(input_data,...
        propn_plot_type,...
        plot_style,...
        colour_vec,...
        yaxis_label,...
        xaxis_label,...
        xticks_vals,...
        xticks_labels,...
        ylim_vals,...
        xlim_vals,...
        display_legend,...
        legend_label,...
        save_filename_prefix,...
        plot_fontsize,....
        dataset)
    
    if ispc==1
        scale = 1/0.75;
    else
        scale = 1;
    end
    
    % Set up plot
    position = [10, 10, 1.5*550*scale, 1.5*450*scale];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    
    % Get number of batches of input data
    n_input_data_batches = numel(input_data);
    
    % Iterate over input data batches
    for inf_itr = 1:n_input_data_batches
        
        % Set violin plotting properties
        if (n_input_data_batches == 2) && (strcmp(dataset,'worker_patterns')==1)
            if inf_itr == 1
                x_offset = -0.15;
            else
                x_offset = 0.15;
            end
            violin_width = 0.1;
        elseif (n_input_data_batches == 2) && (strcmp(dataset,'workpercent_and_backwardsCT')==1)
            if inf_itr == 1
                x_offset = -0.2;
            else
                x_offset = 0.2;
            end
            violin_width = 0.1;
        else
            violin_width = 0.3;
            x_offset = 0;
        end
        
        if strcmp(plot_style,'boxplot')
             % Generate boxplots
             h = boxplot(input_data{inf_itr},...
                            'positions', xticks_vals + x_offset,...
                             'Colors',colour_vec(inf_itr,:),...
                            'Labels',xticks_labels);
             set(h,{'linew'},{1.5})
        elseif strcmp(plot_style,'violin')
            % Generate violin plots
                n_configs = size(input_data{inf_itr},2);
            for config_itr = 1:n_configs

                % If all data points are equal, Violin call will error
                % Plot a line instead
                    config_data = input_data{inf_itr}(:,config_itr);
                unique_vals = unique(config_data);  
                if numel(unique_vals) == 1  % Check if a unique value. If so plot at that point
                    unique_val = unique_vals(1);
                    plot([xticks_vals(config_itr)-0.3 xticks_vals(config_itr)+0.3],...
                          [unique_val unique_val],...
                            'Linewidth',1.5,...
                            'Color',colour_vec(inf_itr,:))
                else
                    if sum(isnan(config_data)) < numel(config_data) % For proportion plots, if baseline values are zero, causes proportional values to all be NaN
                        % Otherwise, create violin plot
                        
                        % Plot violins
                        violins = Violin(config_data, xticks_vals(config_itr) + x_offset,...
                                            'Width',violin_width);

                        % Set violin plot properties
                        violins.ViolinColor = colour_vec(inf_itr,:);
                        violins.ViolinAlpha = 1.0; %shading transparency

                        %Set violin plot region properties
                        violins.EdgeColor = [1 1 1];

                        %Set median marker properties
                        violins.MedianColor = [1 1 1];
                        violins.MedianPlot.Marker = 's';
                        violins.MedianPlot.MarkerEdgeColor = [0 0 0];
                        violins.MedianPlot.LineWidth = 1;

                        %Set whisker line properties
                        violins.BoxColor = [0 0 0];
                    end
                end
            end

            % Set x-axis tick labels
            xticks(xticks_vals)
            xticklabels(xticks_labels)
        else
           error('Incompatible plot type string passed to function.') 
        end
    end
    
    % Set y-axis labels
    ylabel(yaxis_label)
    
    % Set x-axis labels
    xlabel(xaxis_label)
    
    % Set axis limits
    xlim(xlim_vals)
    ylim(ylim_vals)
    
    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on
    
    % Add legend, if applicable
    if display_legend == true
        H = gobjects(n_input_data_batches,1);
        for data_itr = 1:n_input_data_batches
            H(data_itr) = fill(NaN,NaN,colour_vec(data_itr,:),...
                                'DisplayName',legend_label{data_itr});
        end
        leg = legend(H,...
                'LineWidth',1.5,...
                'FontSize',floor(plot_fontsize));
        
            if (strcmp(dataset,'worker_patterns')==1)
                set(leg,'Location','northwest')
            end
    end
    
    
    % Save figure to file
    if strcmp(plot_style,'boxplot')
        if propn_plot_type == true
           save_filename = strcat(save_filename_prefix,'_boxplot_propn');
        else
           save_filename = strcat(save_filename_prefix,'_boxplot');
        end
        %export_fig(save_filename,'-transparent','-pdf','-r1200')
    elseif strcmp(plot_style,'violin')
        if propn_plot_type == true
           save_filename = strcat(save_filename_prefix,'_violins_propn');
        else
           save_filename = strcat(save_filename_prefix,'_violins');
        end
        export_fig(save_filename,'-pdf','-r1200')
    end
end