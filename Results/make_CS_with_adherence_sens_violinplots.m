% PURPOSE:
% Violin plots for COVID-secure workplace outputs.
% Plots relative outputs (to scenario with no COVID-secure workplace
% measures) for different test, trace & isolate adherence assumptions.
% Fixes team size or transmission risk scaling.
%--------------------------------------------------------------------------

function make_CS_with_adherence_sens_violinplots(metric,variable_names)

%% ADD PATH DEPENDENCIES
addpath('Violinplot-Matlab')

%% Load data
% COVID-secure workplace intervention runs.
CS_data = load('worker_model_output_CS_intervention_combined.mat'); % 70% adherence to test, trace & isolate.
CS_data_no_isol = load('worker_model_output_CS_intervention_no_isol_combined.mat'); % 0% adherence to test, trace & isolate.
CS_data_full_isol = load('worker_model_output_CS_intervention_full_isol_combined.mat'); % 0% adherence to test, trace & isolate.

% Also load no CS data (adherence runs)
adherence_data = load('worker_model_output_adherence_intervention_combined.mat');

%% Specify global params

% Specify network size
cmax = 10000;

% Plot fontsize
plot_fontsize = 28;

%% Compute desired summary statistics
if sum(strcmp(variable_names,'final_size') == 1)

    % Set timestep to start counting from
    offset_idx = 15;

    % Number of infections over duration of outbreak
    final_size_absolute_with_isol = squeeze(CS_data.numinf_combined(end,:,:) - CS_data.numinf_combined(offset_idx,:,:));
    final_size_absolute_no_isol = squeeze(CS_data_no_isol.numinf_combined(end,:,:) - CS_data_no_isol.numinf_combined(offset_idx,:,:));
    final_size_absolute_full_isol = squeeze(CS_data_full_isol.numinf_combined(end,:,:) - CS_data_full_isol.numinf_combined(offset_idx,:,:));

    % Propn of infections over duration of outbreak
    input_data_CS_with_isol = final_size_absolute_with_isol/cmax;
    input_data_CS_no_isol = final_size_absolute_no_isol/cmax;
    input_data_CS_full_isol = final_size_absolute_full_isol/cmax;

    % If relative plot, divide by central estimate from no CS runs
    adherence_input_data = squeeze(adherence_data.numinf_combined(end,:,:) - adherence_data.numinf_combined(offset_idx,:,:));
    input_data_no_CS_with_isol = median(adherence_input_data(:,8))/cmax;
    input_data_no_CS_no_isol = median(adherence_input_data(:,1))/cmax;
    input_data_no_CS_full_isol = median(adherence_input_data(:,11))/cmax;

    % Package data for function input
    input_data = {input_data_CS_with_isol,input_data_CS_no_isol,input_data_CS_full_isol,...
                    input_data_no_CS_with_isol,input_data_no_CS_no_isol,input_data_no_CS_full_isol};

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative attack rate';
        ylim_vals = [0 1.5];
    elseif metric == "absolute"
        yaxis_label = 'Additional fraction of the network infected';
        ylim_vals = [0 0.95];
    end

    % Call plot function
    add_legend = true;
    save_filename = 'covid_secure_plots/final_size_CS_adherence_sens_violins';
    variablename = 'final_size';
    plot_grouped_violins_CS_adherence_sense(input_data,...
                                            metric,...
                                            yaxis_label,...
                                            ylim_vals,...
                                            add_legend,...
                                            plot_fontsize,...
                                            save_filename,...
                                            variablename)
end

if sum(strcmp(variable_names,'peak_inf') == 1)

    % Number of infections over duration of outbreak
    input_data_CS_with_isol = squeeze(max(CS_data.prevpresymp_combined +...
                                    CS_data.prevsymp_combined +...
                                    CS_data.prevasymp_combined))/cmax;
    input_data_CS_no_isol = squeeze(max(CS_data_no_isol.prevpresymp_combined +...
                                    CS_data_no_isol.prevsymp_combined +...
                                    CS_data_no_isol.prevasymp_combined))/cmax;
    input_data_CS_full_isol = squeeze(max(CS_data_full_isol.prevpresymp_combined +...
                                    CS_data_full_isol.prevsymp_combined +...
                                    CS_data_full_isol.prevasymp_combined))/cmax;

    % If relative plot, divide by central estimate from no CS runs
    adherence_input_data = squeeze(max(adherence_data.prevpresymp_combined + adherence_data.prevsymp_combined + adherence_data.prevasymp_combined))/cmax;
    input_data_no_CS_with_isol = median(adherence_input_data(:,8));
    input_data_no_CS_no_isol = median(adherence_input_data(:,1));
    input_data_no_CS_full_isol = median(adherence_input_data(:,11));

    % Package data for function input
    input_data = {input_data_CS_with_isol,input_data_CS_no_isol,input_data_CS_full_isol,...
                    input_data_no_CS_with_isol,input_data_no_CS_no_isol,input_data_no_CS_full_isol};

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative peak in infectious case prevalence';
        ylim_vals = [0 1.5];
    elseif metric == "absolute"
        yaxis_label = 'Peak in infectious case prevalence';
        ylim_vals = [0 0.5];
    end

    % Call plot function
    save_filename = 'covid_secure_plots/peak_inf_CS_adherence_sens_violins';
    add_legend = false;
    variablename = 'peak_inf';
    plot_grouped_violins_CS_adherence_sense(input_data,...
                                            metric,...
                                            yaxis_label,...
                                            ylim_vals,...
                                            add_legend,...
                                            plot_fontsize,...
                                            save_filename,...
                                            variablename)
end

if sum(strcmp(variable_names,'duration') == 1)

    % Peak isolations
    prev_data_with_isol = CS_data.prevpresymp_combined + CS_data.prevasymp_combined + CS_data.prevsymp_combined;
    prev_data_no_isol = CS_data_no_isol.prevpresymp_combined + CS_data_no_isol.prevasymp_combined + CS_data_no_isol.prevsymp_combined;
    prev_data_full_isol = CS_data_full_isol.prevpresymp_combined + CS_data_full_isol.prevasymp_combined + CS_data_full_isol.prevsymp_combined;
    input_data_CS_with_isol = zeros(size(squeeze(prev_data_with_isol(1,:,:))));
    input_data_CS_no_isol = zeros(size(squeeze(prev_data_no_isol(1,:,:))));
    input_data_CS_full_isol = zeros(size(squeeze(prev_data_no_isol(1,:,:))));
    for i = 1:length(CS_data.numinf_combined(1,:,1))
        for j = 1:length(CS_data.numinf_combined(1,1,:))
            input_data_CS_with_isol(i,j) = find(prev_data_with_isol(:,i,j)>0,1,'last') - 1;
            input_data_CS_no_isol(i,j) = find(prev_data_no_isol(:,i,j)>0,1,'last') - 1;
            input_data_CS_full_isol(i,j) = find(prev_data_full_isol(:,i,j)>0,1,'last') - 1;
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
    input_data_no_CS_with_isol = median(adherence_duration_data(:,8));
    input_data_no_CS_no_isol = median(adherence_duration_data(:,1));
    input_data_no_CS_full_isol = median(adherence_duration_data(:,11));

    % Package data for function input
    input_data = {input_data_CS_with_isol,input_data_CS_no_isol,input_data_CS_full_isol,...
                    input_data_no_CS_with_isol,input_data_no_CS_no_isol,input_data_no_CS_full_isol};

    % Set y-axis label and limits
    if metric == "relative"
        yaxis_label = 'Relative outbreak duration length';
        ylim_vals = [0 3];
    elseif metric == "absolute"
        yaxis_label = 'Outbreak duration';
        ylim_vals = [0 365];
    end

    % Call plot function
    save_filename = 'covid_secure_plots/durations_CS_adherence_sens_violins';
    add_legend = false;
    variablename = 'duration';
    plot_grouped_violins_CS_adherence_sense(input_data,...
                                            metric,...
                                            yaxis_label,...
                                            ylim_vals,...
                                            add_legend,...
                                            plot_fontsize,...
                                            save_filename,...
                                            variablename)
end

end

%% Function to generate the violin/boxplots plots
function plot_grouped_violins_CS_adherence_sense(input_data,...
                                            metric,...
                                            yaxis_label,...
                                            ylim_vals,...
                                            add_legend,...
                                            plot_fontsize,...
                                            save_filename_prefix,...
                                            variablename)
    % Will create two separate figures
    % First with a fixed team size, varying transmission risk scaling
    % Second with fixed transmission risk scaling, varying team size

    % Unpack input data
    input_data_CS_with_isol = input_data{1};
    input_data_CS_no_isol = input_data{2};
    input_data_CS_full_isol = input_data{3};
    input_data_no_CS_with_isol = input_data{4};
    input_data_no_CS_no_isol = input_data{5};
    input_data_no_CS_full_isol = input_data{6};

    % Amend data based on metric requested
    if metric == "relative"
        input_data_with_isol = input_data_CS_with_isol/input_data_no_CS_with_isol;
        input_data_no_isol = input_data_CS_no_isol/input_data_no_CS_no_isol;
        input_data_full_isol = input_data_CS_full_isol/input_data_no_CS_full_isol;
    elseif metric == "absolute"
        input_data_with_isol = input_data_CS_with_isol;
        input_data_no_isol = input_data_CS_no_isol;
        input_data_full_isol = input_data_CS_full_isol;
    end

    % Set up colour values
    colour_vec = [0.85    0.85   0.85;
                  0.625    0.625   0.625;
                  0.4    0.4   0.4];

    %----------------------------------------------------------------------
    % First figure
    % Fixed team size, varying transmission risk scaling
    position = [10, 10, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(1);
    clf;
    set(fig,'Color', [1 1 1])

    % Get data for a particular fixed work team size
    work_team_size_idxs = 5:8;
    selected_input_data_with_isol = input_data_with_isol(:,work_team_size_idxs);
    selected_input_data_no_isol = input_data_no_isol(:,work_team_size_idxs);
    selected_input_data_full_isol = input_data_full_isol(:,work_team_size_idxs);

    % Set up x-axis ticks
    xticks_vals = 1:2:7;
    xticks_labels = {'0.25','0.50','0.75','1.00'};

    % Set up offset for violin plots based on transmission scaling
    x_offset = [-0.5 0 0.5];

    % Generate violin plots
    % 4 configurations:
    % - Batches of 2 for each CS config
    % - transrisk scaling of 0.25,0.5,0.75,1.0
    for config_itr = 1:4

        % Get relevant data for this iteration
        config_data_with_isol = selected_input_data_with_isol(:,config_itr);
        config_data_no_isol = selected_input_data_no_isol(:,config_itr);
        config_data_full_isol = selected_input_data_full_isol(:,config_itr);

        % Plot violins
        violins_no_isol = Violin(config_data_no_isol, xticks_vals(config_itr) + x_offset(1),...
            'Width',0.2);
        violins_with_isol = Violin(config_data_with_isol, xticks_vals(config_itr) + x_offset(2),...
            'Width',0.2);
        violins_full_isol = Violin(config_data_full_isol, xticks_vals(config_itr) + x_offset(3),...
            'Width',0.2);

        % Set violin plot properties
        violins_no_isol.ViolinColor = colour_vec(1,:);
        violins_no_isol.ViolinAlpha = 1.0; %shading transparency

        violins_with_isol.ViolinColor = colour_vec(2,:);
        violins_with_isol.ViolinAlpha = 1.0; %shading transparency

        violins_full_isol.ViolinColor = colour_vec(3,:);
        violins_full_isol.ViolinAlpha = 1.0; %shading transparency

        %Set violin plot region properties
        violins_with_isol.EdgeColor = [1 1 1];
        violins_no_isol.EdgeColor = [1 1 1];
        violins_full_isol.EdgeColor = [1 1 1];

        %Set median marker properties
        violins_with_isol.MedianColor = [1 1 1];
        violins_with_isol.MedianPlot.Marker = 's';
        violins_with_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_with_isol.MedianPlot.LineWidth = 1;

        violins_no_isol.MedianColor = [1 1 1];
        violins_no_isol.MedianPlot.Marker = 's';
        violins_no_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_no_isol.MedianPlot.LineWidth = 1;

        violins_full_isol.MedianColor = [1 1 1];
        violins_full_isol.MedianPlot.Marker = 's';
        violins_full_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_full_isol.MedianPlot.LineWidth = 1;

        %Set whisker line properties
        violins_with_isol.BoxColor = [0 0 0];
        violins_no_isol.BoxColor = [0 0 0];
        violins_full_isol.BoxColor = [0 0 0];

        % If plotting absolute values, add horizontal line for median value
        % from reference runs with no COVID-secure workplace controls
        if metric == "absolute"
            % 0% adherence
            plot([xticks_vals(config_itr) + x_offset(1) - 0.2 xticks_vals(config_itr) + x_offset(1) + 0.2],...
                    [input_data_no_CS_no_isol input_data_no_CS_no_isol],...
                    '-',...
                    'LineWidth',3.,...
                    'Color',[1 0. 0.]);

            % 70% adherence
            plot([xticks_vals(config_itr) + x_offset(2) - 0.2 xticks_vals(config_itr) + x_offset(2) + 0.2],...
                    [input_data_no_CS_with_isol input_data_no_CS_with_isol],...
                    '-',...
                    'LineWidth',3,...
                    'Color',[1 0. 0.]);

            % 100% adherence
            plot([xticks_vals(config_itr) + x_offset(3) - 0.2 xticks_vals(config_itr) + x_offset(3) + 0.2],...
                    [input_data_no_CS_full_isol input_data_no_CS_full_isol],...
                    '-',...
                    'LineWidth',3,...
                    'Color',[1 0. 0.]);
        end
    end

    % Add horizontal line at median values with
    % no CS measures in use
    if metric == "relative"
        % If relative, add line at y=1
        plot([xticks_vals(1)-1 xticks_vals(end)+1],...
                    [1 1],...
                    '--',...
                    'LineWidth',3.5,...
                    'Color',[0.5 0.5 0.5]);
%     else
%         plot([xticks_vals(1)-1 xticks_vals(end)+1],...
%                     [input_data_no_CS input_data_no_CS],...
%                     '--',...
%                     'LineWidth',3.5,...
%                     'Color',[0.5 0.5 0.5]);
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
    xlabel('Transmission risk scaling')

    % Set axis limits
    xlim([0. 8])
    ylim(ylim_vals)

    % Add legend, if applicable
    legend_label = {'0% adherence','70% adherence','100% adherence'};
    if add_legend == true
        H = gobjects(3,1);
        for data_itr = 1:3
            H(data_itr) = fill(NaN,NaN,colour_vec(data_itr,:),...
                                'DisplayName',legend_label{data_itr});
        end
         leg = legend(H,...
                    'LineWidth',1.5,...
                    'FontSize',26,...
                    'Orientation','vertical',...
                    'Location','northwest');
    end
    %leg.Title.String = 'Transmission risk scaling';
    %leg.Title.FontSize = 20;

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on

    % Save figure to file
    if metric == "absolute"
        leg.Orientation = 'horizontal';
        leg.FontSize = 23;
        save_filename = strcat(save_filename_prefix,'_by_risk');
    else
        save_filename = strcat(save_filename_prefix,'_by_risk_rel');
    end
    export_fig(save_filename,'-pdf','-r1200')

    %----------------------------------------------------------------------
    % Second figure
    % Fixed transmission risk scaling, varying team size
    position = [10, 10, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(2);
    clf;
    set(fig,'Color', [1 1 1])

    % Get data for a particular fixed transmission risk
    risk_idxs = [2,6,10]; % [2,6,10], transmission risk of 50%
    selected_input_data_with_isol = input_data_with_isol(:,risk_idxs);
    selected_input_data_no_isol = input_data_no_isol(:,risk_idxs);
    selected_input_data_full_isol = input_data_full_isol(:,risk_idxs);

    % Set up x-axis ticks
    xticks_vals = 1:2:5;
    xticks_labels = {'2','5','10'};

    % Set up offset for violin plots based on transmission scaling
    x_offset = [-0.5 0 0.5];

    % Generate violin plots
    % 3 configurations:
    % - Batches of 2 for each CS config
    % - work team size of 2, 5, 10
    for config_itr = 1:3

         % Get relevant data for this iteration
        config_data_with_isol = selected_input_data_with_isol(:,config_itr);
        config_data_no_isol = selected_input_data_no_isol(:,config_itr);
        config_data_full_isol = selected_input_data_full_isol(:,config_itr);

        % Plot violins
        violins_no_isol = Violin(config_data_no_isol, xticks_vals(config_itr) + x_offset(1),...
            'Width',0.2);
        violins_with_isol = Violin(config_data_with_isol, xticks_vals(config_itr) + x_offset(2),...
            'Width',0.2);
        violins_full_isol = Violin(config_data_full_isol, xticks_vals(config_itr) + x_offset(3),...
            'Width',0.2);

        % Set violin plot properties
        violins_no_isol.ViolinColor = colour_vec(1,:);
        violins_no_isol.ViolinAlpha = 1.0; %shading transparency

        violins_with_isol.ViolinColor = colour_vec(2,:);
        violins_with_isol.ViolinAlpha = 1.0; %shading transparency

        violins_full_isol.ViolinColor = colour_vec(3,:);
        violins_full_isol.ViolinAlpha = 1.0; %shading transparency

        %Set violin plot region properties
        violins_with_isol.EdgeColor = [1 1 1];
        violins_no_isol.EdgeColor = [1 1 1];
        violins_full_isol.EdgeColor = [1 1 1];

        %Set median marker properties
        violins_with_isol.MedianColor = [1 1 1];
        violins_with_isol.MedianPlot.Marker = 's';
        violins_with_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_with_isol.MedianPlot.LineWidth = 1;

        violins_no_isol.MedianColor = [1 1 1];
        violins_no_isol.MedianPlot.Marker = 's';
        violins_no_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_no_isol.MedianPlot.LineWidth = 1;

        violins_full_isol.MedianColor = [1 1 1];
        violins_full_isol.MedianPlot.Marker = 's';
        violins_full_isol.MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_full_isol.MedianPlot.LineWidth = 1;

        %Set whisker line properties
        violins_with_isol.BoxColor = [0 0 0];
        violins_no_isol.BoxColor = [0 0 0];
        violins_full_isol.BoxColor = [0 0 0];

        % If plotting absolute values, add horizontal line for median value
        % from reference runs with no COVID-secure workplace controls
        if metric == "absolute"
            % 0% adherence
            plot([xticks_vals(config_itr) + x_offset(1) - 0.2 xticks_vals(config_itr) + x_offset(1) + 0.2],...
                    [input_data_no_CS_no_isol input_data_no_CS_no_isol],...
                    '-',...
                    'LineWidth',3,...
                    'Color',[1 0. 0.]);

            % 70% adherence
            plot([xticks_vals(config_itr) + x_offset(2) - 0.2 xticks_vals(config_itr) + x_offset(2) + 0.2],...
                    [input_data_no_CS_with_isol input_data_no_CS_with_isol],...
                    '-',...
                    'LineWidth',3,...
                    'Color',[1 0. 0.]);

            % 100% adherence
            plot([xticks_vals(config_itr) + x_offset(3) - 0.2 xticks_vals(config_itr) + x_offset(3) + 0.2],...
                    [input_data_no_CS_full_isol input_data_no_CS_full_isol],...
                    '-',...
                    'LineWidth',3,...
                    'Color',[1 0. 0.]);
        end
    end

    % Add horizontal line at median values with
    % no CS measures in use
    if metric == "relative"
        % If relative, add line at y=1
        plot([xticks_vals(1)-1 xticks_vals(end)+1],...
                    [1 1],...
                    '--',...
                    'LineWidth',3.5,...
                    'Color',[0.5 0.5 0.5]);
%     else
%         plot([xticks_vals(1)-1 xticks_vals(end)+1],...
%                     [input_data_no_CS input_data_no_CS],...
%                     '--',...
%                     'LineWidth',3.5,...
%                     'Color',[0.5 0.5 0.5]);
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
    xlim([0 6])
    ylim(ylim_vals)

%     % Add legend, if applicable
%     legend_label = {'With other interventions','No other interventions'};
%     if add_legend == true
%         H = gobjects(2,1);
%         for data_itr = 1:2
%             H(data_itr) = fill(NaN,NaN,colour_vec(data_itr,:),...
%                                 'DisplayName',legend_label{data_itr});
%         end
%         leg = legend(H,...
%                 'LineWidth',1.5,...
%                 'FontSize',22,...
%                 'Orientation','vertical',...
%                 'Location','northwest');
%     end
%     %leg.Title.String = 'Transmission risk scaling';
%     %leg.Title.FontSize = 20;

    %Specify general axis properties
    set(gca,'FontSize',floor(plot_fontsize))
    set(gca,'LineWidth',1)
    box on

    % Save figure to file
    if metric == "absolute"
        save_filename = strcat(save_filename_prefix,'_by_team_size');
    else
        save_filename = strcat(save_filename_prefix,'_by_team_size_rel');
    end
    export_fig(save_filename,'-pdf','-r1200')
    %----------------------------------------------------------------------
end
