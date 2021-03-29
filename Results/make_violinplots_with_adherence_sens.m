% PURPOSE:
% Make violin plots for 0% and 100% adherence runs
%--------------------------------------------------------------------------

% clear variables

function make_violinplots_with_adherence_sens(dataset,variablename)

%% ADD PATH DEPENDENCIES
addpath('Violinplot-Matlab')

%% Load data
% Work percent
workpercent_data_no_isol = load('worker_model_output_workpercent_intervention_low_adherence_combined.mat');
workpercent_data_full_isol = load('worker_model_output_workpercent_intervention_full_adherence_combined.mat');

% Synchronised
synch_data_no_isol = load('worker_model_output_synchronised_changedays_intervention_low_adherence_combined.mat');
synch_data_full_isol = load('worker_model_output_synchronised_changedays_intervention_full_adherence_combined.mat');

% Non-synchronised
asynch_data_no_isol = load('worker_model_output_variable_changedays_intervention_low_adherence_combined.mat');
asynch_data_full_isol = load('worker_model_output_variable_changedays_intervention_full_adherence_combined.mat');

%% Specify global params

% Specify network size
cmax = 10000;

%% Compute desired summary statistics
if strcmp(variablename,'final_size')
    % Number of infections from day 15 onwards (i.e. after day 14, which is timestep 15)
    offset_idx = 15;
    workpercent_no_isol_final_size = squeeze(workpercent_data_no_isol.numinf_combined(end,:,:) - workpercent_data_no_isol.numinf_combined(offset_idx,:,:));
    workpercent_full_isol_final_size = squeeze(workpercent_data_full_isol.numinf_combined(end,:,:) - workpercent_data_full_isol.numinf_combined(offset_idx,:,:));

    synch_no_isol_final_size = squeeze(synch_data_no_isol.numinf_combined(end,:,:) - synch_data_no_isol.numinf_combined(offset_idx,:,:));
    synch_full_isol_final_size = squeeze(synch_data_full_isol.numinf_combined(end,:,:) - synch_data_full_isol.numinf_combined(offset_idx,:,:));

    asynch_no_isol_final_size = squeeze(asynch_data_no_isol.numinf_combined(end,:,:) - asynch_data_no_isol.numinf_combined(offset_idx,:,:));
    asynch_full_isol_final_size = squeeze(asynch_data_full_isol.numinf_combined(end,:,:) - asynch_data_full_isol.numinf_combined(offset_idx,:,:));

    % Proportion of infections from day 15 onwards
    workpercent_no_isol_input_data = workpercent_no_isol_final_size/cmax;
    workpercent_full_isol_input_data = workpercent_full_isol_final_size/cmax;

    synch_no_isol_input_data = synch_no_isol_final_size/cmax;
    synch_full_isol_input_data = synch_full_isol_final_size/cmax;

    asynch_no_isol_input_data = asynch_no_isol_final_size/cmax;
    asynch_full_isol_input_data = asynch_full_isol_final_size/cmax;
elseif strcmp(variablename,'peak_inf')
    % Peak proportion infections
    workpercent_no_isol_input_data = squeeze(max(workpercent_data_no_isol.prevpresymp_combined(:,:,:) + workpercent_data_no_isol.prevsymp_combined(:,:,:) + workpercent_data_no_isol.prevasymp_combined(:,:,:)))/cmax;
    workpercent_full_isol_input_data = squeeze(max(workpercent_data_full_isol.prevpresymp_combined(:,:,:) + workpercent_data_full_isol.prevsymp_combined(:,:,:) + workpercent_data_full_isol.prevasymp_combined(:,:,:)))/cmax;

    synch_no_isol_input_data = squeeze(max(synch_data_no_isol.prevpresymp_combined + synch_data_no_isol.prevsymp_combined + synch_data_no_isol.prevasymp_combined))/cmax;
    synch_full_isol_input_data = squeeze(max(synch_data_full_isol.prevpresymp_combined + synch_data_full_isol.prevsymp_combined + synch_data_full_isol.prevasymp_combined))/cmax;

    asynch_no_isol_input_data = squeeze(max(asynch_data_no_isol.prevpresymp_combined + asynch_data_no_isol.prevsymp_combined + asynch_data_no_isol.prevasymp_combined))/cmax;
    asynch_full_isol_input_data = squeeze(max(asynch_data_full_isol.prevpresymp_combined + asynch_data_full_isol.prevsymp_combined + asynch_data_full_isol.prevasymp_combined))/cmax;
elseif strcmp(variablename,'total_isolation')
    % total number of isolation-days from day 15 onwards
    offset_idx = 15;
    workpercent_no_isol_input_data = squeeze(sum(workpercent_data_no_isol.num_isolating_combined(offset_idx:end,:,:)));
    workpercent_full_isol_input_data = squeeze(sum(workpercent_data_full_isol.num_isolating_combined(offset_idx:end,:,:)));

    synch_no_isol_input_data = squeeze(sum(synch_data_no_isol.num_isolating_combined(offset_idx:end,:,:)));
    synch_full_isol_input_data = squeeze(sum(synch_data_full_isol.num_isolating_combined(offset_idx:end,:,:)));

    asynch_no_isol_input_data = squeeze(sum(asynch_data_no_isol.num_isolating_combined(offset_idx:end,:,:)));
    asynch_full_isol_input_data = squeeze(sum(asynch_data_full_isol.num_isolating_combined(offset_idx:end,:,:)));
elseif strcmp(variablename,'duration')
    % outbreak duration
    workpercent_data_no_isol.prev = workpercent_data_no_isol.prevpresymp_combined(:,:,:) + workpercent_data_no_isol.prevasymp_combined(:,:,:) + workpercent_data_no_isol.prevsymp_combined(:,:,:);
    workpercent_data_full_isol.prev = workpercent_data_full_isol.prevpresymp_combined(:,:,:) + workpercent_data_full_isol.prevasymp_combined(:,:,:) + workpercent_data_full_isol.prevsymp_combined(:,:,:);

    synch_data_no_isol.prev = synch_data_no_isol.prevpresymp_combined + synch_data_no_isol.prevasymp_combined + synch_data_no_isol.prevsymp_combined;
    synch_data_full_isol.prev = synch_data_full_isol.prevpresymp_combined + synch_data_full_isol.prevasymp_combined + synch_data_full_isol.prevsymp_combined;

    asynch_data_no_isol.prev = asynch_data_no_isol.prevpresymp_combined + asynch_data_no_isol.prevasymp_combined + asynch_data_no_isol.prevsymp_combined;
    asynch_data_full_isol.prev = asynch_data_full_isol.prevpresymp_combined + asynch_data_full_isol.prevasymp_combined + asynch_data_full_isol.prevsymp_combined;

    workpercent_no_isol_input_data = zeros(size(squeeze(workpercent_data_no_isol.prev(1,:,:))));
    workpercent_full_isol_input_data = zeros(size(squeeze(workpercent_data_full_isol.prev(1,:,:))));

    synch_no_isol_input_data = zeros(size(squeeze(synch_data_no_isol.prev(1,:,:))));
    synch_full_isol_input_data = zeros(size(squeeze(synch_data_full_isol.prev(1,:,:))));

    asynch_no_isol_input_data = zeros(size(squeeze(asynch_data_no_isol.prev(1,:,:))));
    asynch_full_isol_input_data = zeros(size(squeeze(asynch_data_full_isol.prev(1,:,:))));

    for i = 1:length(workpercent_data_no_isol.numinf_combined(1,:,1))
        for j = 1:length(workpercent_data_no_isol.numinf_combined(1,1,:))
            workpercent_no_isol_input_data(i,j) = find(workpercent_data_no_isol.prev(:,i,j)>0,1,'last') - 1;
            workpercent_full_isol_input_data(i,j) = find(workpercent_data_full_isol.prev(:,i,j)>0,1,'last') - 1;
        end
    end
    for i = 1:length(synch_data_no_isol.numinf_combined(1,:,1))
        for j = 1:length(synch_data_no_isol.numinf_combined(1,1,:))
            synch_no_isol_input_data(i,j) = find(synch_data_no_isol.prev(:,i,j)>0,1,'last') - 1;
            synch_full_isol_input_data(i,j) = find(synch_data_full_isol.prev(:,i,j)>0,1,'last') - 1;

            asynch_no_isol_input_data(i,j) = find(asynch_data_no_isol.prev(:,i,j)>0,1,'last') - 1;
            asynch_full_isol_input_data(i,j) = find(asynch_data_full_isol.prev(:,i,j)>0,1,'last') - 1;
        end
    end
end

%% WOrkpercent sensitivity
if strcmp(dataset,'workpercent')==1

    % Specify input data
    input_data = {workpercent_no_isol_input_data,workpercent_full_isol_input_data};

    % Set plot style
    plot_style = 'violin';

    % Set colour for violin plots
    colour_vec = [0.85 0.85 0.85;
                    0.5 0.5 0.5];

    % Set x-axis label
    xaxis_label = 'Proportion working from home';

    % Set up xticks
    xticks_vals = 0:1:11;
    xticks_labels = {'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','N-U'};

    % Set up xaxis limits
    xlim_vals = [-0.5 11.5];

    % Set filename prefix
    save_filename_initial = ['workpercent_plots/worker_model_sens_workpercent_and_adherence_',variablename];

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
    legend_label = {'10% adherence','100% adherence'};
elseif strcmp(dataset,'synch')==1
    %% Synch & asynch worker pattern comparison

    % Specify input data
    input_data = {synch_no_isol_input_data,synch_full_isol_input_data};

    % Set plot style
    plot_style = 'violin';

    % Set colour for violin plots
    colour_vec = [0.85 0.85 0.85;
                    0.5 0.5 0.5];

    % Set x-axis label
    xaxis_label = 'Days per week at the workplace';

    % Set up xticks
    xticks_vals = 0:1:5;
    xticks_labels = {'0','1','2','3','4','5'};

    % Set up xaxis limits
    xlim_vals = [-0.5 5.5];

    % Set filename prefix
    save_filename_initial = ['worker_pattern_plots/worker_model_sens_synch_worker_pattern_and_adherence_',variablename];

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
    legend_label = {'10% adherence','100% adherence'};
elseif strcmp(dataset,'asynch')==1
    %% Synch & asynch worker pattern comparison

    % Specify input data
    input_data = {asynch_no_isol_input_data,asynch_full_isol_input_data};

    % Set plot style
    plot_style = 'violin';

    % Set colour for violin plots
    colour_vec = [0.85 0.85 0.85;
                    0.5 0.5 0.5];

    % Set x-axis label
    xaxis_label = 'Days per week at the workplace';

    % Set up xticks
    xticks_vals = 0:1:5;
    xticks_labels = {'0','1','2','3','4','5'};

    % Set up xaxis limits
    xlim_vals = [-0.5 5.5];

    % Set filename prefix
    save_filename_initial = ['worker_pattern_plots/worker_model_sens_asynch_worker_pattern_and_adherence_',variablename];

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
    legend_label = {'10% adherence','100% adherence'};
end

% Set up y-axis labels per statistic.
% If required, normalise the data
if propn_plot_type == true
    if strcmp(variablename,'final_size')
        yaxis_label = 'Relative proportion infected';
    elseif strcmp(variablename,'peak_inf')
        yaxis_label = 'Relative peak in infectious case prevalence';
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
        yaxis_label = 'Peak in infectious case prevalence';
    elseif strcmp(variablename,'total_isolation')
        yaxis_label = 'Total isolation-days';
    elseif strcmp(variablename,'peak_isolation')
        yaxis_label = 'Peak proportion in isolation';
    elseif strcmp(variablename,'duration')
        yaxis_label = 'Outbreak duration';
    end
end

% Set up yaxis limits
if strcmp(variablename,'total_isolation')
    ylim_vals = [0 350000];
elseif strcmp(variablename,'final_size')
    ylim_vals = [0 0.9];
elseif strcmp(variablename,'peak_inf')
    ylim_vals = [0 0.55];
elseif strcmp(variablename,'duration')
    ylim_vals = [0 370];
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
        dataset,...
        variablename)
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
        dataset,...
        variablename)


    % Set up plot
    position = [10, 10, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])

    % Get number of batches of input data
    n_input_data_batches = numel(input_data);

    % Iterate over input data batches
    for inf_itr = 1:n_input_data_batches

        % Set violin plotting properties
        if (n_input_data_batches == 2)
            if inf_itr == 1
                x_offset = -0.15;
            else
                x_offset = 0.15;
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
    if strcmp(variablename,'final_size')
        ytickformat('%.1f')
    elseif (strcmp(variablename,'peak_inf'))
        ytickformat('%.2f')
    elseif (strcmp(variablename,'total_isolation'))
        ytickformat('%.1f')
    end

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

            if (strcmp(dataset,'synch')==1)
                set(leg,'Location','northwest')
            elseif (strcmp(dataset,'asynch')==1)
                set(leg,'Location','northwest')
            elseif (strcmp(dataset,'workpercent')==1)
                set(leg,'Location','northeast','Orientation','horizontal')
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
