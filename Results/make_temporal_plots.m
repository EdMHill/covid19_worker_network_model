function maxy = make_temporal_plots(dataset,variablename,make_subplot)
% does a temporal plot of the variablename, for each dataset in the list of datasets
% if make_subplot==1 then a new figure isn't generated, and the plot isn't
% saved
% e.g. make_temporal_plots({'adherence','workpercent'},'inf_prevalence',0)


%% Load data
for i=1:length(dataset)
    if strcmp(dataset{i},'adherence')
        % Adherence
        data{i} = load('worker_model_output_adherence_intervention_combined.mat');
    elseif strcmp(dataset{i},'workpercent')
        % Work percent
        data{i} = load('worker_model_output_workpercent_intervention_combined.mat');
    elseif strcmp(dataset{i},'backwards_CT')
        % Backwards CT
        data{i} = load('worker_model_output_amount_backwards_CT_intervention_combined.mat');
    elseif strcmp(dataset{i},'synch')
        % Synchronised worker pattern
        data{i} = load('worker_model_output_synchronised_changedays_intervention_combined.mat');
    elseif strcmp(dataset{i},'asynch')
        % Asynchronised worker pattern
        data{i} = load('worker_model_output_variable_changedays_intervention_combined.mat');
    end
end

%% Specify global params

% Specify network size
cmax = 10000;

%% Compute desired summary statistics
for i=1:length(data)
    if strcmp(variablename,'inf_prevalence')
        % Infectious prevalence
        temporal_data{i} = (data{i}.prevpresymp_combined + data{i}.prevsymp_combined + data{i}.prevasymp_combined)/cmax;
    elseif strcmp(variablename,'isol_prevalence')
        % Isolation prevalence
        temporal_data{i} = (data{i}.num_isolating_combined)/cmax;
    elseif strcmp(variablename,'Rt')
        % Rt window size
        Rt_window_size = 7;
       
        % Get data
        Rt_data = data{i}.Rt_save_combined;
        temporal_data{i} = movmean(Rt_data,Rt_window_size,'omitnan','Endpoints','discard');
    end
    
    % Get quantiles at each timestep across the profiles
    % Entry breakdown: 1 - median, 2&3 - 50% PI, 4&5 - 90% PI, 6&7 - 95%,
    % 8&9 - 99%
    prctile_vals{i} = [50 25 75 5 95 2.5 97.5 0.5 99.5];
    temporal_data_prctiles{i} = prctile(temporal_data{i},prctile_vals{i},2);
    temporal_data_median{i} = squeeze(temporal_data_prctiles{i}(:,1,:));
    temporal_data_95PI_lower{i} = squeeze(temporal_data_prctiles{i}(:,6,:));
    temporal_data_95PI_upper{i} = squeeze(temporal_data_prctiles{i}(:,7,:));
end

% Set up y-axis labels per statistic.
% If required, normalise the data
if strcmp(variablename,'inf_prevalence')
    yaxis_label = 'Propn. infectious';
elseif strcmp(variablename,'isol_prevalence')
    yaxis_label = 'Propn. isolating';
end

% Set up xticks
x_vals = 0:365;
x_cutoff = 16;
if strcmp(variablename,'Rt')
    if mod(Rt_window_size,2) == 0
        x_vals_Rt = (floor(Rt_window_size/2)):(365 - floor(Rt_window_size/2)+1);
    else
        x_vals_Rt = (floor(Rt_window_size/2)):(365 - floor(Rt_window_size/2));
    end
    x_cutoff_Rt = x_cutoff - floor(Rt_window_size/2);
end

% Set up xaxis limits
xlim_vals = [0 365];
maxy = max(max(temporal_data_prctiles{1}(:,3,:)));
for i=1:length(data)
    maxy = max([maxy, max(max(temporal_data_prctiles{1}(:,3,:)))]);
end
ylim_vals = [0 maxy*1.05];

%% Plot temporal profile of number of infectious prevalence
% Have panel per team size and transmission risk

% Set up plot colours
grey_colour_vec = [0.5 0.5 0.5];
% Set colour for temporal plots for each dataset
for i=1:length(data)
    if strcmp(dataset{i},'adherence')
        % Adherence
        % create a default color map ranging from base colour to much
        % lighter
        num_lines{i} = length(temporal_data_prctiles{i}(1,1,:));
        basecolour = [0.8500    0.3250    0.0980];
        white = [1 1 1];
        colors_p = [linspace(basecolour(1),white(1),num_lines{i}+3)', linspace(basecolour(2),white(2),num_lines{i}+3)', linspace(basecolour(3),white(3),num_lines{i}+3)'];
        colour_vec{i} = colors_p(1:num_lines{i},:);
    elseif strcmp(dataset{i},'workpercent')
        % Work percent
        % create a default color map ranging from base colour to much
        % lighter
        num_lines{i} = 11;
        basecolour = [0.6 0 0.8];
        white = [1 1 1];
        colors_p = [linspace(basecolour(1),white(1),num_lines{i}+2)', linspace(basecolour(2),white(2),num_lines{i}+2)', linspace(basecolour(3),white(3),num_lines{i}+2)'];
        colour_vec{i} = colors_p(1:num_lines{i},:);
    elseif strcmp(dataset{i},'backwards_CT')
        % Backwards contact tracing
        % create a default color map ranging from base colour to much
        % lighter
        num_lines{i} = 8;
        basecolour = [0 0 1];
        white = [1 1 1];
        colors_p = [linspace(basecolour(1),white(1),num_lines{i}+3)', linspace(basecolour(2),white(2),num_lines{i}+3)', linspace(basecolour(3),white(3),num_lines{i}+3)'];
        colour_vec{i} = colors_p(1:num_lines{i},:);
    elseif strcmp(dataset{i},'synch')
        % Synchronised worker pattern
        % create a default color map ranging from base colour to much
        % lighter
        num_lines{i} = 6;
        basecolour = [0.4 0 0];
        white = [1 1 1];
        colors_p = [linspace(white(1),basecolour(1),num_lines{i}+2)', linspace(white(2),basecolour(2),num_lines{i}+2)', linspace(white(3),basecolour(3),num_lines{i}+2)'];
        colour_vec{i} = colors_p(1:num_lines{i},:);
    elseif strcmp(dataset{i},'asynch')
        % Asynchronised worker pattern
        % create a default color map ranging from base colour to much
        % lighter
        num_lines{i} = 6;
        basecolour = [0.2 0.8 0.8];
        white = [1 1 1];
        colors_p = [linspace(white(1),basecolour(1),num_lines{i}+2)', linspace(white(2),basecolour(2),num_lines{i}+2)', linspace(white(3),basecolour(3),num_lines{i}+2)'];
        colour_vec{i} = colors_p(1:num_lines{i},:);    
    end
end

if make_subplot~=1
    % Set up plot
    if ispc==1
        scale = 1/0.75;
    else
        scale = 1;
    end
    position = [10, 10, 3.5*550*scale, 4.5*450*scale];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    hold on
end

% Construct plots. For Rt need to use different x-values to plot against
if strcmp(variablename,'Rt')
%     % Fill the 50% prediction intervals up to day x_cutoff
%     fill([x_vals_Rt(1:end) x_vals_Rt(x_cutoff_Rt:-1:1)], [temporal_data_prctiles{i}(1:x_cutoff_Rt,2,1); temporal_data_prctiles{i}(x_cutoff_Rt:-1:1,3,1)],'r','FaceColor',1-0.7*(1-grey_colour_vec),'EdgeColor','none');
% 
%     % Plot traces up to day x_cutoff
%     plot(x_vals_Rt(1:x_cutoff_Rt),temporal_data_median{i}(1:x_cutoff_Rt,1),'Linewidth',1.5,'Color',grey_colour_vec)

    for i=1:length(data)

        for array_idx = num_lines{i}:-1:1

%             % Fill the 50% prediction intervals after day x_cutoff
%             fill([x_vals_Rt(1:end) x_vals_Rt(end:-1:1)], [temporal_data_prctiles{i}(1:end,2,array_idx); temporal_data_prctiles{i}(end:-1:1,3,array_idx)],'r','FaceColor',1-0.7*(1-colour_vec{i}(array_idx,:)),'EdgeColor','none');

            % Plot traces after day x_cutoff
            plot(x_vals_Rt(1:end),temporal_data_median{i}(1:end,array_idx),'Linewidth',1.5,'Color',colour_vec{i}(array_idx,:))
%             
%             % Fill the 50% prediction intervals after day x_cutoff
%             fill([x_vals_Rt(x_cutoff_Rt:end) x_vals_Rt(end:-1:x_cutoff_Rt)], [temporal_data_prctiles{i}(x_cutoff_Rt:end,2,array_idx); temporal_data_prctiles{i}(end:-1:x_cutoff_Rt,3,array_idx)],'r','FaceColor',1-0.7*(1-colour_vec{i}(array_idx,:)),'EdgeColor','none');
% 
%             % Plot traces after day x_cutoff
%             plot(x_vals_Rt(x_cutoff_Rt:end),temporal_data_median{i}(x_cutoff_Rt:end,array_idx),'Linewidth',1.5,'Color',colour_vec{i}(array_idx,:))
        end
    end
else

    % Fill the 50% prediction intervals up to day x_cutoff
    fill([x_vals(1:x_cutoff) x_vals(x_cutoff:-1:1)], [temporal_data_prctiles{i}(1:x_cutoff,2,1); temporal_data_prctiles{i}(x_cutoff:-1:1,3,1)],'r','FaceColor',1-0.7*(1-grey_colour_vec),'EdgeColor','none');

    % Plot traces up to day x_cutoff
    plot(x_vals(1:x_cutoff),temporal_data_median{i}(1:x_cutoff,1),'Linewidth',1.5,'Color',grey_colour_vec)

    for i=1:length(data)

        for array_idx = num_lines{i}:-1:1
            % Fill the 50% prediction intervals after day x_cutoff
            fill([x_vals(x_cutoff:end) x_vals(end:-1:x_cutoff)], [temporal_data_prctiles{i}(x_cutoff:end,2,array_idx); temporal_data_prctiles{i}(end:-1:x_cutoff,3,array_idx)],'r','FaceColor',1-0.7*(1-colour_vec{i}(array_idx,:)),'EdgeColor','none');

            % Plot traces after day x_cutoff
            plot(x_vals(x_cutoff:end),temporal_data_median{i}(x_cutoff:end,array_idx),'Linewidth',1.5,'Color',colour_vec{i}(array_idx,:))
        end
    end
end

% set transparency so you can see the lines (alpha between 0 (transparent
% and 1 opaque)
alpha 0.25

%Specify axes properties
xlim(xlim_vals)
ylim(ylim_vals)

%Specify general axis properties
set(gca,'FontSize',26)
set(gca,'LineWidth',1)
box on


% % Add legend to top left panel
% subplot(plot_nRows,plot_nCols,plot_nCols)
% p1 = plot(NaN,NaN,'-','Color',1-0.7*(1-colour_vec),'LineWidth',1.5,'DisplayName','Median');
% p2 = fill(NaN,NaN,'r','FaceColor',1-0.7*(1-colour_vec),'EdgeColor','none','DisplayName','50% prediction interval');
% p3 = fill(NaN,NaN,'r','FaceColor',1-0.45*(1-colour_vec),'EdgeColor','none','DisplayName','90% prediction interval');
% p4 = fill(NaN,NaN,'r','FaceColor',1-0.2*(1-colour_vec),'EdgeColor','none','DisplayName','99% prediction interval');
% legend([p1;p2;p3;p4],...
%     'FontSize',20,...
%     'Position',[0.79 0.86 0.093 0.04]);

if make_subplot~=1
    ylabel(yaxis_label)
    xlabel('Time (days)')
    % Save file
    save_filename = ['misc_plots/temporal_profiles_',variablename];
    for i=1:length(data)
        save_filename = [save_filename,'_',dataset{i}];
    end
    export_fig(save_filename,'-pdf','-r600')
end