% PURPOSE:
% Violin plots for COVID-secure workplace outputs
% Group of three per transmission scaling. Violin plot per team size
%--------------------------------------------------------------------------

function make_CS_violinplots(dataset,metric,variable_names)

%% ADD PATH DEPENDENCIES
addpath('Violinplot-Matlab')

%% Load data
if strcmp(dataset,'COVID_secure')==1
    CS_data = load('worker_model_output_CS_intervention_combined.mat');
elseif strcmp(dataset,'COVID_secure_no_isol')==1
    CS_data = load('worker_model_output_CS_intervention_no_isol_combined.mat');
end

% Also load no CS data (adherence runs)
adherence_data = load('worker_model_output_adherence_intervention_combined.mat');

%% Specify global params

% Specify network size
cmax = 10000;

% Number of simulations replicates performed
n_simns = 1000;

% Plot fontsize
plot_fontsize = 24;

%% Specify save file identifier
if strcmp(dataset,'COVID_secure')==1
    save_file_addon = '_with_isol';
elseif strcmp(dataset,'COVID_secure_no_isol')==1
    save_file_addon = '_no_isol';
end

%% Compute desired summary statistics
if sum(strcmp(variable_names,'final_size') == 1)

    % Set timestep to start counting from
    offset_idx = 15;

    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(CS_data.numinf_combined(end,:,:) - CS_data.numinf_combined(offset_idx,:,:));

    % Propn of infections over duration of outbreak
    input_data = final_size_absolute/cmax;

    % If relative plot, divide by central estimate from no CS runs
    adherence_input_data = squeeze(adherence_data.numinf_combined(end,:,:) - adherence_data.numinf_combined(offset_idx,:,:));
    if strcmp(dataset,'COVID_secure')==1
        input_data_no_CS = median(adherence_input_data(:,8))/cmax;
    elseif strcmp(dataset,'COVID_secure_no_isol')==1
        input_data_no_CS = median(adherence_input_data(:,1))/cmax;
    end

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative attack rate';
        ylim_vals = [0 1.5];
    elseif metric == "absolute"
        yaxis_label = 'Additional fraction of the network infected';
        ylim_vals = [0 0.85];
    end

    % Call plot function
    if strcmp(dataset,'COVID_secure')==1
        add_legend = true;
    else
        add_legend = false;
    end
    save_filename = ['covid_secure_plots/final_size_violins',save_file_addon];
    variablename = 'final_size';
    plot_grouped_violins(input_data,input_data_no_CS,metric,...
                            yaxis_label,...
                            ylim_vals,...
                            add_legend,plot_fontsize,save_filename,...
                            variablename)

%     % Check against the threshold criteria and store
%     heatmap_vals = sum(final_size_propn>0.5)/n_simns;
%
%     % Reorder into a 2D array. Row by transrisk, column by worker group size
%     data_array = reshape(heatmap_vals,4,3);
%
%     % Call plot function
%     yaxis_label = 'Additional fraction of the network infected';
%     add_legend = true;
%     save_filename = ['covid_secure_plots/final_size_violins',save_file_addon];
%     plot_grouped_violins(heatmap_array,yaxis_label,xaxis_label,panel_title,...
%                             add_colourbar,plot_fontsize,save_filename)
end

if sum(strcmp(variable_names,'peak_inf') == 1)

    % Number of infections over duration of outbreak
    input_data = squeeze(max(CS_data.prevpresymp_combined +...
                                    CS_data.prevsymp_combined +...
                                    CS_data.prevasymp_combined))/cmax;

    % If relative plot, divide by central estimate from no CS runs
    adherence_input_data = squeeze(max(adherence_data.prevpresymp_combined + adherence_data.prevsymp_combined + adherence_data.prevasymp_combined))/cmax;
    if strcmp(dataset,'COVID_secure')==1
        input_data_no_CS = median(adherence_input_data(:,8));
    elseif strcmp(dataset,'COVID_secure_no_isol')==1
        input_data_no_CS = median(adherence_input_data(:,1));
    end

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative peak in infectious case prevalence';
        ylim_vals = [0 1.5];
    elseif metric == "absolute"
        yaxis_label = 'Peak in infectious case prevalence';
        ylim_vals = [0 0.3];
    end

    % Call plot function
    save_filename = ['covid_secure_plots/peak_inf_violins',save_file_addon];
    add_legend = false;
    variablename = 'peak_inf';
    plot_grouped_violins(input_data,input_data_no_CS,metric,...
                            yaxis_label,...
                            ylim_vals,...
                            add_legend,plot_fontsize,save_filename,...
                            variablename)
end

if sum(strcmp(variable_names,'total_isolation') == 1)

    % Time in isolation (average per person)
    offset_idx = 15;
    input_data = squeeze(sum(CS_data.num_isolating_combined(offset_idx:end,:,:)));

    % If relative plot, divide by central estimate from no CS runs
    adherence_input_data = squeeze(sum(adherence_data.num_isolating_combined(offset_idx:end,:,:)));
    if strcmp(dataset,'COVID_secure')==1
        input_data_no_CS = median(adherence_input_data(:,8));
    elseif strcmp(dataset,'COVID_secure_no_isol')==1
        input_data_no_CS = median(adherence_input_data(:,1));
    end

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative amount of isolation-days';
        ylim_vals = [0 1.5];
    elseif metric == "absolute"
        yaxis_label = 'Total isolation-days';
        ylim_vals = [0 19e4];
    end

    % Call plot function
    save_filename = ['covid_secure_plots/avg_isol_violins',save_file_addon];
    add_legend = false;
    variablename = 'total_isolation';
    plot_grouped_violins(input_data,input_data_no_CS,metric,...
                            yaxis_label,...
                            ylim_vals,...
                            add_legend,plot_fontsize,save_filename,...
                            variablename)
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

    % If relative plot, divide by central estimate from no CS runs
    adherence_prev_data = adherence_data.prevpresymp_combined + adherence_data.prevasymp_combined + adherence_data.prevsymp_combined;
    adherence_duration_data = zeros(size(squeeze(adherence_prev_data(1,:,:))));
    for i = 1:length(adherence_data.numinf_combined(1,:,1))
        for j = 1:length(adherence_data.numinf_combined(1,1,:))
            adherence_duration_data(i,j) = find(adherence_prev_data(:,i,j)>0,1,'last') - 1;
        end
    end

    if strcmp(dataset,'COVID_secure')==1
        input_data_no_CS = median(adherence_duration_data(:,8));
    elseif strcmp(dataset,'COVID_secure_no_isol')==1
        input_data_no_CS = median(adherence_duration_data(:,1));
    end

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative outbreak duration length';
        ylim_vals = [0 3];
    elseif metric == "absolute"
        yaxis_label = 'Outbreak duration';
        ylim_vals = [0 365];
    end

    % Call plot function
    save_filename = ['covid_secure_plots/durations_violins',save_file_addon];
    add_legend = false;
    variablename = 'duration';
    plot_grouped_violins(duration_data,input_data_no_CS,metric,...
                            yaxis_label,...
                            ylim_vals,...
                            add_legend,plot_fontsize,save_filename,...
                            variablename)
end

end

%% Function to generate the violin/boxplots plots
function  plot_grouped_violins(input_data,...
                                input_data_no_CS,...
                                metric,...
                                yaxis_label,...
                                ylim_vals,...
                                add_legend,...
                                plot_fontsize,...
                                save_filename_prefix,...
                                variablename)

    % Amend data based on metric requested
    if metric == "relative"
        input_data = input_data/input_data_no_CS;
    end

    % Set up colour values
    colour_vec = [0    0.4470    0.7410;
                  0.8500    0.3250    0.0980;
                  0.9290    0.6940    0.1250;
                  0.75    0.4    0.8];

%     % Set up plot
%     position = [10, 10, 1.5*550, 1.5*450];
%     set(0, 'DefaultFigurePosition', position);
%     %set(0, 'DefaultFigurePosition', position,'Units','points');
%     fig = figure(1);
%     clf;
%     set(fig,'Color', [1 1 1])
%
%
%     % Set up x-axis ticks
%     xticks_vals = 0:1:3;
%     xticks_labels = {'0.25','0.50','0.75','1.0'};
%
%     % Set up offset for violin plots based on worker team size
%     x_offset = [-0.3 0 0.3];
%
%     % Generate violin plots
%     % 12 configurations:
%     % - Batches of 4 for each worker team size
%     % - transrisk scaling of 0.25,0.5,0.75,1.0
%     n_configs = size(input_data,2);
%     for config_itr = 1:n_configs
%
%         % Get worker size batch being analysed
%         worker_size_itr = ceil(config_itr/4);
%
%         % Get transmission scaling being analysed
%         transmission_scale_itr = mod(config_itr,4);
%         if transmission_scale_itr == 0
%             transmission_scale_itr = 4;
%         end
%
%         % If all data points are equal, Violin call will error
%         % Plot a line instead
%         config_data = input_data(:,config_itr);
%         unique_vals = unique(config_data);
%         if numel(unique_vals) == 1  % Check if a unique value. If so plot at that point
%             unique_val = unique_vals(1);
%             plot([xticks_vals(transmission_scale_itr)-0.3 xticks_vals(transmission_scale_itr)+0.3],...
%                 [unique_val unique_val],...
%                 'Linewidth',1.5,...
%                 'Color',colour_vec(worker_size_itr,:))
%         else
%             if sum(isnan(config_data)) < numel(config_data) % For proportion plots, if baseline values are zero, causes proportional values to all be NaN
%                 % Otherwise, create violin plot
%
%                 % Plot violins
%                 violins = Violin(config_data, xticks_vals(transmission_scale_itr) + x_offset(worker_size_itr),...
%                     'Width',0.15);
%
%                 % Set violin plot properties
%                 violins.ViolinColor = colour_vec(worker_size_itr,:);
%                 violins.ViolinAlpha = 1.0; %shading transparency
%
%                 %Set violin plot region properties
%                 violins.EdgeColor = [1 1 1];
%
%                 %Set median marker properties
%                 violins.MedianColor = [1 1 1];
%                 violins.MedianPlot.Marker = 's';
%                 violins.MedianPlot.MarkerEdgeColor = [0 0 0];
%                 violins.MedianPlot.LineWidth = 1;
%
%                 %Set whisker line properties
%                 violins.BoxColor = [0 0 0];
%             end
%         end
%     end
%
%     % If using absolute counts, add horizontal line at median values with
%     % no CS measures in use
%     if metric == "absolute"
%         plot([xticks_vals(1)-0.5 xticks_vals(end)+0.5],...
%                     [input_data_no_CS input_data_no_CS],...
%                     '--',...
%                     'LineWidth',1.5,...
%                     'Color',[0.5 0.5 0.5]);
%     end
%
%     % Set x-axis tick labels
%     xticks(xticks_vals)
%     xticklabels(xticks_labels)
%
%     % Set y-axis labels
%     ylabel(yaxis_label)
%     ytickformat('%.2f')
%
%     % Set x-axis labels
%     xlabel('Work contact transmission risk scaling')
%
% %     % Set axis limits
% %     xlim(xlim_vals)
% %     ylim(ylim_vals)
%
%     %Specify general axis properties
%     set(gca,'FontSize',floor(plot_fontsize))
%     set(gca,'LineWidth',1)
%     box on
% %
% %     % Add legend, if applicable
% %     if display_legend == true
% %         H = gobjects(n_input_data_batches,1);
% %         for data_itr = 1:n_input_data_batches
% %             H(data_itr) = fill(NaN,NaN,colour_vec(data_itr,:),...
% %                                 'DisplayName',legend_label{data_itr});
% %         end
% %         leg = legend(H,...
% %                 'LineWidth',1.5,...
% %                 'FontSize',floor(plot_fontsize));
% %     end
%
%     % Save figure to file
%     save_filename = strcat(save_filename_prefix,'_violins');
%     %export_fig(save_filename,'-pdf','-r1200')

    % Set up plot
    position = [10, 10, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure(2);
    clf;
    set(fig,'Color', [1 1 1])

    % Set up x-axis ticks
    xticks_vals = 0:1.5:3;
    xticks_labels = {'2','5','10'};

    % Set up offset for violin plots based on transmission scaling
    x_offset = [-0.3 -0.1 0.1 0.3];

    % Generate violin plots
    % 12 configurations:
    % - Batches of 4 for each worker team size
    % - transrisk scaling of 0.25,0.5,0.75,1.0
    n_configs = size(input_data,2);
    for config_itr = 1:n_configs

        % Get worker size batch being analysed
        worker_size_itr = ceil(config_itr/4);

        % Get transmission scaling being analysed
        transmission_scale_itr = mod(config_itr,4);
        if transmission_scale_itr == 0
            transmission_scale_itr = 4;
        end

        % If all data points are equal, Violin call will error
        % Plot a line instead
        config_data = input_data(:,config_itr);
        if sum(isnan(config_data)) < numel(config_data) % For proportion plots, if baseline values are zero, causes proportional values to all be NaN
            % Plot violins
            violins = Violin(config_data, xticks_vals(worker_size_itr) + x_offset(transmission_scale_itr),...
                'Width',0.075);

            % Set violin plot properties
            violins.ViolinColor = colour_vec(transmission_scale_itr,:);
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

    % If using absolute counts, add horizontal line at median values with
    % no CS measures in use
    if metric == "absolute"
        plot([xticks_vals(1)-0.5 xticks_vals(end)+0.5],...
                    [input_data_no_CS input_data_no_CS],...
                    '--',...
                    'LineWidth',3.5,...
                    'Color',[0.5 0.5 0.5]);
    else
        % If relative, add line at y=1
        plot([xticks_vals(1)-0.5 xticks_vals(end)+0.5],...
                    [1 1],...
                    '--',...
                    'LineWidth',3.5,...
                    'Color',[0.5 0.5 0.5]);
    end

    % Set x-axis tick labels
    xticks(xticks_vals)
    xticklabels(xticks_labels)

    % Set y-axis labels
    ylabel(yaxis_label)
    if metric == "relative"
        ytickformat('%.1f')
    else
        if strcmp(variablename,'final_size')
            ytickformat('%.1f')
        elseif (strcmp(variablename,'peak_inf'))
            ytickformat('%.2f')
        elseif (strcmp(variablename,'total_isolation'))
            ytickformat('%.0f')
        end
    end

    % Set x-axis labels
    xlabel('Work team size')

    % Set axis limits
    xlim([-0.45 3.45])
    ylim(ylim_vals)

    % Add legend, if applicable
    legend_label = {'0.25','0.50','0.75','1.00'};
    if add_legend == true
        H = gobjects(4,1);
        for data_itr = 1:4
            H(data_itr) = fill(NaN,NaN,colour_vec(data_itr,:),...
                                'DisplayName',legend_label{data_itr});
        end
        leg = legend(H,...
                'LineWidth',1.5,...
                'FontSize',20,...
                'Orientation','horizontal',...
                'Location','northwest');
    end
    leg.Title.String = 'Transmission risk scaling';
    leg.Title.FontSize = 20;

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on

    % Save figure to file
    if metric == "absolute"
        save_filename = strcat(save_filename_prefix,'_violins_by_team_size');
    else
        save_filename = strcat(save_filename_prefix,'_violins_by_team_size_rel');
    end
    export_fig(save_filename,'-pdf','-r1200')
end
