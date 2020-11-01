% PURPOSE:
% Produce figures based on the synchronised and unsynchronised work
% pattern runs
%--------------------------------------------------------------------------

clear variables

%% Load data
% Synchronised
synch_data = load('../worker_model_output_synchronised_changedays_intervention_combined.mat');

% Non-synchronised
unsynch_data = load('../worker_model_output_variable_changedays_intervention_combined.mat');

%% Compute desired summary statistics

% Specify network size
cmax = 10000;

% Infectious prevalence
synch_data_infectious_prev = (synch_data.prevpresymp_combined + synch_data.prevsymp_combined + synch_data.prevasymp_combined)/cmax;
unsynch_data_infectious_prev = (unsynch_data.prevpresymp_combined + unsynch_data.prevsymp_combined + unsynch_data.prevasymp_combined)/cmax;

    % Get quantiles at each timestep across the profiles
    % Entry breakdown: 1 - median, 2&3 - 50% PI, 4&5 - 90% PI, 6&7 - 95%,
    % 8&9 - 99%
prctile_vals = [50 25 75 5 95 2.5 97.5 0.5 99.5];   
synch_data_infectious_prev_prctiles = prctile(synch_data_infectious_prev,prctile_vals,2);
unsynch_data_infectious_prev_prctiles = prctile(unsynch_data_infectious_prev,prctile_vals,2);

synch_data_infectious_prev_median = squeeze(synch_data_infectious_prev_prctiles(:,1,:));
unsynch_data_infectious_prev_median = squeeze(unsynch_data_infectious_prev_prctiles(:,1,:));

%% Plot temporal profile of number of infectious prevalence
%% First: synchronised work patterns; Second: Variable work patterns

% Set global plot properties
y_max = 400/cmax;
x_max = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First: synchronised work patterns %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up plot
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
%set(0, 'DefaultFigurePosition', position,'Units','points');
fig = figure();
clf;
set(fig,'Color', [1 1 1])
hold on

% Plot traces
plot(synch_data_infectious_prev_median,'Linewidth',1.5)
%plot(squeeze(synch_data_infectious_prev_prctiles(31:end,6,:)),'--','Linewidth',1.5)
%plot(squeeze(synch_data_infectious_prev_prctiles(31:end,7,:)),'--','Linewidth',1.5)

%Specify axes properties
xlim([0 x_max])
ylim([0 y_max])

% Set axes labels
xlabel('Day')
ylabel('Proportion infectious')

% Add plot title
title('Synchronised')

%Specify general axis properties
set(gca,'FontSize',26)
set(gca,'LineWidth',1)
box on

% Save file
export_fig('worker_pattern_sync_temporal_medians','-pdf','-transparent','-r1200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second: Variable work patterns %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up plot
position = [10, 10, 1.5*550, 1.5*450];
set(0, 'DefaultFigurePosition', position);
%set(0, 'DefaultFigurePosition', position,'Units','points');
fig = figure();
clf;
set(fig,'Color', [1 1 1])
hold on

% Plot traces
plot(unsynch_data_infectious_prev_median,'Linewidth',1.5)

%Specify axes properties
xlim([0 x_max])
ylim([0 y_max])

% Set axes labels
xlabel('Day')
ylabel('Proportion infectious')

% Add plot title
title('Asynchronised')

% Add legend
legend({'Full lockdown','1 day','2 days','3 days','4 days','5 days'})

%Specify general axis properties
set(gca,'FontSize',26)
set(gca,'LineWidth',1)
box on

% Save file
export_fig('worker_pattern_async_temporal_medians','-pdf','-transparent','-r1200')