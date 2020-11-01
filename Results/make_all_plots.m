
%% make final size plots
make_violinplots('adherence','final_size')
make_violinplots('workpercent_and_backwardsCT','final_size')

%% make peak infections plots
make_violinplots('adherence','peak_inf')
make_violinplots('workpercent_and_backwardsCT','peak_inf')

%% total self-isolating days
make_violinplots('adherence','total_isolation')
make_violinplots('workpercent_and_backwardsCT','total_isolation')

%% peak fraction self-isolating
make_violinplots('adherence','peak_isolation')
make_violinplots('workpercent_and_backwardsCT','peak_isolation')

%% outbreak duration
make_violinplots('adherence','duration')
make_violinplots('workpercent_and_backwardsCT','duration')

%% Threshold event plots: Adherence
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation','duration'};
criteria_thresholds = [0.2,0.01,0.01,0.05,150];
display_legend = true;
make_threshold_event_plots('adherence',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Backwards CT
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation'};
criteria_thresholds = [0.2,0.01,0.01,0.05];
display_legend = false;
make_threshold_event_plots('backwardsCT',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Workpercent
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation'};
criteria_thresholds = [0.2,0.01,0.01,0.05];
display_legend = true;
make_threshold_event_plots('workpercent',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Backwards CT & Workpercent
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation','duration'};
criteria_thresholds = [0.2,0.01,0.01,0.05,150];
display_legend = true;
make_threshold_event_plots('workpercent_and_backwardsCT',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Synch & asynch worker patterns
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation','duration'};
criteria_thresholds = [0.2,0.01,0.01,0.05,150];
display_legend = true;
make_threshold_event_plots('worker_patterns',variable_names,criteria_thresholds,display_legend)

%% COVID-secure temporal outputs
make_temporal_grid('COVID_secure','inf_prevalence')
make_temporal_grid('COVID_secure','isol_prevalence')
make_temporal_grid('COVID_secure_no_isol','inf_prevalence')
%% COVID-secure heatmap outputs
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation','duration'};
make_heatmaps('COVID_secure',variable_names)
make_heatmaps('COVID_secure_no_isol',variable_names)