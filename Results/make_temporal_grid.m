% PURPOSE:
% Produce grid to show variability in outputs across range of values for
% two input variables.
%--------------------------------------------------------------------------

function make_temporal_grid(dataset,variablename)


%% Load data
if strcmp(dataset,'COVID_secure')==1
    CS_data = load('worker_model_output_CS_intervention_combined.mat');
elseif strcmp(dataset,'COVID_secure_no_isol')==1
    CS_data = load('worker_model_output_CS_intervention_no_isol_combined.mat');
end
%% Specify global params

% Specify network size
cmax = 10000;

%% Compute desired summary statistics
if strcmp(variablename,'inf_prevalence')
    % Infectious prevalence
    CS_temporal_data = (CS_data.prevpresymp_combined + CS_data.prevsymp_combined + CS_data.prevasymp_combined)/cmax;
elseif strcmp(variablename,'isol_prevalence')
    % Isolation prevalence
    CS_temporal_data = CS_data.num_isolating_combined/(cmax*365); 
end

% Get quantiles at each timestep across the profiles
% Entry breakdown: 1 - median, 2&3 - 50% PI, 4&5 - 90% PI, 6&7 - 95%,
% 8&9 - 99%
prctile_vals = [50 25 75 5 95 2.5 97.5 0.5 99.5];   
CS_temporal_data_prctiles = prctile(CS_temporal_data,prctile_vals,2);
CS_temporal_data_median = squeeze(CS_temporal_data_prctiles(:,1,:));
CS_temporal_data_99PI_upper = squeeze(CS_temporal_data_prctiles(:,9,:));

%% Specify plotting variabels based on dataset
if (strcmp(dataset,'COVID_secure')==1) || (strcmp(dataset,'COVID_secure_no_isol')==1)
    
    % Set up xticks
    x_vals = 0:365;
    x_cutoff = 30;
    
    % Set axis limits
    xlim_vals = [0 200];
    ylim_vals = [0 max(CS_temporal_data_99PI_upper(:))*1.05];
    
    % Set plot fontsize
    plot_fontsize = 26;
    
    % Set variable values
    panel_xVals = [2,5,10];          % Team size
    panel_yVals = [0.25,0.5,0.75,1]; % Transrisk scaling
    
    % Set number of columns and rows in plots
    plot_nCols = numel(panel_xVals);
    plot_nRows = numel(panel_yVals);
end

% Set up y-axis labels per statistic.
% If required, normalise the data
if strcmp(variablename,'inf_prevalence')
    yaxis_label = 'Proportion infectious';
    colour_vec = [0.8500    0.3250    0.0980];
elseif strcmp(variablename,'isol_prevalence')
    yaxis_label = 'Proportion isolating';
    colour_vec = [0    0    1];
end

%% Plot temporal profile of number of infectious prevalence
% Have panel per team size and transmission risk

% Set up plot colours
grey_colour_vec = [0.5 0.5 0.5];

% Set y-axis values
if (strcmp(variablename,'inf_prevalence')) && (strcmp(dataset,'COVID_secure'))
    ytick_vals = [0,0.05,0.10,0.15];
elseif (strcmp(variablename,'isol_prevalence')) && (strcmp(dataset,'COVID_secure'))
    ytick_vals = [0 0.0002 0.0004  0.0006 0.0008 0.001];
elseif (strcmp(variablename,'inf_prevalence')) && (strcmp(dataset,'COVID_secure_no_isol'))
    ytick_vals = [0,0.05,0.10,0.15];
end

% Set up plot
position = [10, 10, 3.5*550, 4.5*450];
set(0, 'DefaultFigurePosition', position);
%set(0, 'DefaultFigurePosition', position,'Units','points');
fig = figure();
clf;
set(fig,'Color', [1 1 1])
hold on

% Initialise panel ID
panel_idx = 1;

for row_itr = 1:plot_nRows
    for col_itr = 1:plot_nCols
        % Specify panel to be plotted to
        subplot(plot_nRows,plot_nCols,panel_idx)
        hold on
        
        % Get relevant scenario from the dataset
        array_idx = row_itr + ((col_itr-1)*plot_nRows);

        % Fill the prediction intervals up to day x_cutoff
        % 99%. 90% & 50%
        fill([x_vals(1:x_cutoff) x_vals(x_cutoff:-1:1)], [CS_temporal_data_prctiles(1:x_cutoff,8,array_idx); CS_temporal_data_prctiles(x_cutoff:-1:1,9,array_idx)],'r','FaceColor',1-0.2*(1-grey_colour_vec),'EdgeColor','none');
        fill([x_vals(1:x_cutoff) x_vals(x_cutoff:-1:1)], [CS_temporal_data_prctiles(1:x_cutoff,4,array_idx); CS_temporal_data_prctiles(x_cutoff:-1:1,5,array_idx)],'r','FaceColor',1-0.45*(1-grey_colour_vec),'EdgeColor','none');
        fill([x_vals(1:x_cutoff) x_vals(x_cutoff:-1:1)], [CS_temporal_data_prctiles(1:x_cutoff,2,array_idx); CS_temporal_data_prctiles(x_cutoff:-1:1,3,array_idx)],'r','FaceColor',1-0.7*(1-grey_colour_vec),'EdgeColor','none');

        % Plot traces up to day x_cutoff
        plot(x_vals(1:x_cutoff),CS_temporal_data_median(1:x_cutoff,array_idx),'Linewidth',1.5,'Color',grey_colour_vec)

        % Fill the prediction intervals after day x_cutoff
        % 99%. 90% & 50%
        fill([x_vals(x_cutoff:end) x_vals(end:-1:x_cutoff)], [CS_temporal_data_prctiles(x_cutoff:end,8,array_idx); CS_temporal_data_prctiles(end:-1:x_cutoff,9,array_idx)],'r','FaceColor',1-0.2*(1-colour_vec),'EdgeColor','none');
        fill([x_vals(x_cutoff:end) x_vals(end:-1:x_cutoff)], [CS_temporal_data_prctiles(x_cutoff:end,4,array_idx); CS_temporal_data_prctiles(end:-1:x_cutoff,5,array_idx)],'r','FaceColor',1-0.45*(1-colour_vec),'EdgeColor','none');
        fill([x_vals(x_cutoff:end) x_vals(end:-1:x_cutoff)], [CS_temporal_data_prctiles(x_cutoff:end,2,array_idx); CS_temporal_data_prctiles(end:-1:x_cutoff,3,array_idx)],'r','FaceColor',1-0.7*(1-colour_vec),'EdgeColor','none');

        % Plot traces after day x_cutoff
        plot(x_vals(x_cutoff:end),CS_temporal_data_median(x_cutoff:end,array_idx),'Linewidth',1.5,'Color',colour_vec)

        %Specify axes properties
        xlim(xlim_vals)
        ylim(ylim_vals)
        
         % Set y-axis ticks
         yticks(ytick_vals)

        % Set axes labels
%         if row_itr == plot_nRows
%             xlabel('Time (Days)')
%         end
%         if col_itr == 1
%             ylabel(yaxis_label)
%         end
%         if strcmp(variablename,'inf_prevalence')
%             if panel_idx == 1
%                 text(-43,-0.315,yaxis_label,'FontSize',plot_fontsize,'FontWeight','bold','Rotation',90)
%             elseif panel_idx == 11
%                 text(60,-0.05,'Time (Days)','FontSize',plot_fontsize,'FontWeight','bold')
%             end
%         end

        %Specify general axis properties
        set(gca,'FontSize',plot_fontsize)
        set(gca,'LineWidth',1)
        box on
        
        % Increment the panel to be plotted to
        panel_idx = panel_idx + 1;
    end
end

% Set headers and yticks
if (strcmp(variablename,'inf_prevalence')) && (strcmp(dataset,'COVID_secure'))
 
    % Set up column headers
    text(-480,0.82,'Work team size: 2','FontSize',plot_fontsize,'FontWeight','bold')
    text(-218,0.82,'Work team size: 5','FontSize',plot_fontsize,'FontWeight','bold')
    text(40,0.82,'Work team size: 10','FontSize',plot_fontsize,'FontWeight','bold')

    % Set up row headers
    text(268,0.65,'COVID-secure transmission risk scaling','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.75,'0.25','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.53,'0.50','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.32,'0.75','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.10,'1.00','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)

    % Set up axis labels
    text(-566,0.3,yaxis_label,'FontSize',plot_fontsize,'FontWeight','bold','Rotation',90)
    text(-193,-0.05,'Time (Days)','FontSize',plot_fontsize,'FontWeight','bold')
elseif (strcmp(variablename,'isol_prevalence')) && (strcmp(dataset,'COVID_secure'))      
    
    % Set up column headers
    text(-480,0.0062,'Work team size: 2','FontSize',plot_fontsize,'FontWeight','bold')
    text(-218,0.0062,'Work team size: 5','FontSize',plot_fontsize,'FontWeight','bold')
    text(40,0.0062,'Work team size: 10','FontSize',plot_fontsize,'FontWeight','bold')

    % Set up row headers
    text(268,0.0048,'COVID-secure transmission risk scaling','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.0057,'0.25','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.0041,'0.50','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.0024,'0.75','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.00068,'1.00','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)

    % Set up axis labels
    text(-560,0.00205,yaxis_label,'FontSize',plot_fontsize,'FontWeight','bold','Rotation',90)
    text(-193,-0.0004,'Time (Days)','FontSize',plot_fontsize,'FontWeight','bold')
elseif (strcmp(variablename,'inf_prevalence')) && (strcmp(dataset,'COVID_secure_no_isol'))             
    
    % Set up column headers
    text(-480,0.95,'Work team size: 2','FontSize',plot_fontsize,'FontWeight','bold')
    text(-218,0.95,'Work team size: 5','FontSize',plot_fontsize,'FontWeight','bold')
    text(40,0.95,'Work team size: 10','FontSize',plot_fontsize,'FontWeight','bold')

    % Set up row headers
    text(268,0.75,'COVID-secure transmission risk scaling','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.86,'0.25','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.62,'0.50','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.37,'0.75','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)
    text(228,0.12,'1.00','FontSize',plot_fontsize,'FontWeight','bold','Rotation',270)

    % Set up axis labels
    text(-566,0.28,yaxis_label,'FontSize',plot_fontsize,'FontWeight','bold','Rotation',90)
    text(-193,-0.05,'Time (Days)','FontSize',plot_fontsize,'FontWeight','bold')
end

% Add legend to top left panel
subplot(plot_nRows,plot_nCols,plot_nCols)
p1 = plot(NaN,NaN,'-','Color',1-0.7*(1-colour_vec),'LineWidth',1.5,'DisplayName','Median');
p2 = fill(NaN,NaN,'r','FaceColor',1-0.7*(1-colour_vec),'EdgeColor','none','DisplayName','50% prediction interval');
p3 = fill(NaN,NaN,'r','FaceColor',1-0.45*(1-colour_vec),'EdgeColor','none','DisplayName','90% prediction interval');
p4 = fill(NaN,NaN,'r','FaceColor',1-0.2*(1-colour_vec),'EdgeColor','none','DisplayName','99% prediction interval');
legend([p1;p2;p3;p4],...
    'FontSize',20,...
    'Position',[0.78 0.835 0.09 0.04]);

% Save file
if strcmp(dataset,'COVID_secure')
    save_filename = ['covid_secure_plots/temporal_profiles_',variablename];
elseif strcmp(dataset,'COVID_secure_no_isol')
    save_filename = ['covid_secure_plots/temporal_profiles_no_isol_',variablename];
end
export_fig(save_filename,'-pdf','-r600')

end

% %% Compute desired summary statistics
% if strcmp(variablename,'final_size')
%     % Number of infections from day 30 onwards (i.e. timestep 31)
%     offset_idx = 31;
%     CS_final_size = squeeze(CS_data.numinf_combined(end,:,:) - CS_data.numinf_combined(offset_idx,:,:));
% 
%     % Number of infections from day 30 onwards (i.e. timestep 31)
%     CS_input_data = CS_final_size/cmax;
% elseif strcmp(variablename,'peak_inf')
%     % Peak proportion new infections
%     CS_input_data = squeeze(max(CS_data.newinf_combined))/cmax;
% elseif strcmp(variablename,'total_isolation')
%     % total number of isolation-days from day 30 onwards (i.e. timestep 31)
%     offset_idx = 31;
%     CS_input_data = squeeze(sum(CS_data.num_isolating_combined(offset_idx:end,:,:)));
% elseif strcmp(variablename,'peak_isolation')
%     % peak number isolating
%     CS_input_data = squeeze(max(CS_data.num_isolating_combined))/cmax;
% end
% 
% % Set up y-axis labels per statistic.
% % If required, normalise the data
% if propn_plot_type == true
%     if strcmp(variablename,'final_size')
%         yaxis_label = 'Relative proportion infected';
%     elseif strcmp(variablename,'peak_inf')
%         yaxis_label = 'Relative peak proportion new infections';
%     elseif strcmp(variablename,'total_isolation')    
%         yaxis_label = 'Relative total isolation-days';    
%     elseif strcmp(variablename,'peak_isolation')
%         yaxis_label = 'Relative peak proportion in isolation';
%     elseif strcmp(variablename,'duration')
%         yaxis_label = 'Relative outbreak duration';
%     end
% else
%     if strcmp(variablename,'final_size')
%         yaxis_label = 'Additional fraction of the network infected';
%     elseif strcmp(variablename,'peak_inf')
%         yaxis_label = 'Peak proportion new infections';
%     elseif strcmp(variablename,'total_isolation')
%         yaxis_label = 'Total isolation-days';  
%     elseif strcmp(variablename,'peak_isolation')
%         yaxis_label = 'Peak proportion in isolation';
%     elseif strcmp(variablename,'duration')
%         yaxis_label = 'Outbreak duration';
%     end
% end
% 
% %% Call plot function
% 
% % If required, normalise values relative to first column of data
% if propn_plot_type == true
%     baseline_vals = input_data(:,1);
%     input_data = input_data./baseline_vals;
% end
% 
% generate_sensitivity_plot(input_data,...
%         propn_plot_type,...
%         plot_style,...
%         colour_vec,...
%         yaxis_label,...
%         xaxis_label,...
%         xticks_vals,...
%         xticks_labels,...
%         ylim_vals,...
%         xlim_vals,...
%         display_legend,...
%         legend_label,...
%         save_filename_initial,...
%         plot_fontsize)
