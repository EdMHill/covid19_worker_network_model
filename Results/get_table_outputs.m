% PURPOSE:
% Get outputs for tables
%--------------------------------------------------------------------------

clear variables


%% Specify global params

% Specify network size
cmax = 10000;

% Number of simulations replicates performed
n_simns = 1000;

%% Load data
adherence_data = load('worker_model_output_adherence_intervention_combined.mat');
workpercent_data = load('worker_model_output_workpercent_intervention_combined.mat');
backwards_CT_data = load('worker_model_output_amount_backwards_CT_intervention_combined.mat');
synch_data = load('worker_model_output_synchronised_changedays_intervention_combined.mat');
asynch_data = load('worker_model_output_variable_changedays_intervention_combined.mat');
CS_data = load('worker_model_output_CS_intervention_combined.mat');
CS_data_no_isol = load('worker_model_output_CS_intervention_no_isol_combined.mat');

%% Threshold event plots
format bank
adherence_array = threshold_event_vals(adherence_data,cmax,n_simns)
workpercent_array = threshold_event_vals(workpercent_data,cmax,n_simns)
backwards_CT_array = threshold_event_vals(backwards_CT_data,cmax,n_simns)
synch_array = threshold_event_vals(synch_data,cmax,n_simns)
asynch_array = threshold_event_vals(asynch_data,cmax,n_simns)

%% Violin plots
format long
adherence_violins_array = violin_vals(adherence_data,cmax,n_simns)
workpercent_violins_array = violin_vals(workpercent_data,cmax,n_simns)
backwards_CT_violins_array = violin_vals(backwards_CT_data,cmax,n_simns)
synch_violins_array = violin_vals(synch_data,cmax,n_simns)
asynch_violins_array = violin_vals(asynch_data,cmax,n_simns)

%% Heatmaps
format bank
CS_data_no_isol_final_size = heatmap_table_vals(CS_data_no_isol,'final_size',cmax,n_simns)
CS_data_no_isol_duration = heatmap_table_vals(CS_data_no_isol,'duration',cmax,n_simns)
CS_data_final_size = heatmap_table_vals(CS_data,'final_size',cmax,n_simns)
CS_data_duration = heatmap_table_vals(CS_data,'duration',cmax,n_simns)
CS_data_isol = heatmap_table_vals(CS_data,'avg_isolation',cmax,n_simns)
CS_data_peak_isol = heatmap_table_vals(CS_data,'peak_isolation',cmax,n_simns)

%% Functions to compute desired summary statistics

function output_array = threshold_event_vals(input_data,cmax,n_simns)

    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(input_data.numinf_combined(end,:,:));

    % Propn of infections over duration of outbreak
    final_size_propn = final_size_absolute/cmax;

    % Check against the threshold criteria and store
    final_size_vals = sum(final_size_propn>0.5)/n_simns;

    % Number of infections over duration of outbreak
    peak_inf = squeeze(max(input_data.newinf_combined))/cmax;

    % Check against the threshold criteria and store
    peak_inf_vals = sum(peak_inf>0.01)/n_simns;

    % Time in isolation (average per person)
    total_isol = squeeze(sum(input_data.num_isolating_combined(1:end,:,:)))/(cmax*365);

    % Check against the threshold criteria and store
    total_isol_vals = sum(total_isol>0.01)/n_simns;

    % Peak isolations
    peak_isol = squeeze(max(input_data.num_isolating_combined))/cmax;

    % Check against the threshold criteria and store
    peak_isol_vals = sum(peak_isol>0.05)/n_simns;

    % Duration
    prev_data = input_data.prevpresymp_combined + input_data.prevasymp_combined + input_data.prevsymp_combined;
    duration_data = zeros(size(squeeze(prev_data(1,:,:))));
    for i = 1:length(input_data.numinf_combined(1,:,1))
        for j = 1:length(input_data.numinf_combined(1,1,:))
            duration_data(i,j) = find(prev_data(:,i,j)>0,1,'last');        
        end
    end
    
    % Check against the threshold criteria and store
    duration_vals = sum(duration_data>150)/n_simns;


    % Put into output array
    output_array = [final_size_vals;peak_inf_vals;total_isol_vals;duration_vals;peak_isol_vals];
end

function output_array = violin_vals(input_data,cmax,n_simns)

    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(input_data.numinf_combined(end,:,:));

    % Propn of infections over duration of outbreak
    final_size_propn = final_size_absolute/cmax;

    % Number of infections over duration of outbreak
    peak_inf = squeeze(max(input_data.newinf_combined))/cmax;

    % Time in isolation (average per person)
    total_isol = squeeze(sum(input_data.num_isolating_combined(1:end,:,:)));

    % Peak isolations
    peak_isol = squeeze(max(input_data.num_isolating_combined))/cmax;

    % Duration
    prev_data = input_data.prevpresymp_combined + input_data.prevasymp_combined + input_data.prevsymp_combined;
    duration_data = zeros(size(squeeze(prev_data(1,:,:))));
    for i = 1:length(input_data.numinf_combined(1,:,1))
        for j = 1:length(input_data.numinf_combined(1,1,:))
            duration_data(i,j) = find(prev_data(:,i,j)>0,1,'last');        
        end
    end
    
    % Set percentile values
    % Entry breakdown: 1 - median, 2&3 - 50% PI, 4&5 - 90% PI, 6&7 - 95%,
    % 8&9 - 99%
    prctile_vals = [50 25 75 5 95 2.5 97.5 0.5 99.5];
    
    % Initialise output array
    n_stats = 5;
    n_scens = size(input_data.numinf_combined,3);
    output_array = cell(n_stats,n_scens);
    for stat_itr = 1:n_stats
        if stat_itr == 1
            summ_stat_data = final_size_propn;
            formatSpec = '%.2f';
        elseif stat_itr == 2
            summ_stat_data = peak_inf;
            formatSpec = '%.2f';
        elseif stat_itr == 3
            summ_stat_data = total_isol;
            formatSpec = '%.0u';
        elseif stat_itr == 4
            summ_stat_data = duration_data;
            formatSpec = '%.0u';
        elseif stat_itr == 5
            summ_stat_data = peak_isol;
            formatSpec = '%.2f';
        end
        
        % Get percentiles for current summary statistic
        data_prctiles = prctile(summ_stat_data,prctile_vals,1);
        
        if stat_itr == 3
            data_prctiles = round(data_prctiles,-2);
        elseif stat_itr == 4
            data_prctiles = round(data_prctiles);
        end

        for scen_itr = 1:n_scens
            output_array{stat_itr,scen_itr} = [num2str(data_prctiles(1,scen_itr),formatSpec),...
                                                ' (',num2str(data_prctiles(6,scen_itr),formatSpec),',',num2str(data_prctiles(7,scen_itr),formatSpec),')'];
        end
    end
end

function heatmap_array = heatmap_table_vals(CS_data,variable_name,cmax,n_simns)
if strcmp(variable_name,'final_size')

    % Number of infections over duration of outbreak
    final_size_absolute = squeeze(CS_data.numinf_combined(end,:,:));

    % Propn of infections over duration of outbreak
    final_size_propn = final_size_absolute/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(final_size_propn>0.5)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
end

if sum(strcmp(variable_name,'peak_inf') == 1)

    % Number of infections over duration of outbreak
    peak_inf_CS_data = squeeze(max(CS_data.newinf_combined))/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(peak_inf_CS_data>0.01)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
end

if sum(strcmp(variable_name,'avg_isolation') == 1)
        
    % Time in isolation (average per person)
    offset_idx = 31;
    total_isol_CS_data = squeeze(sum(CS_data.num_isolating_combined(offset_idx:end,:,:)))/(cmax*365);

    % Check against the threshold criteria and store
    heatmap_vals = sum(total_isol_CS_data>0.01)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
end

if sum(strcmp(variable_name,'peak_isolation') == 1)

    % Peak isolations
    peak_isol_CS_data = squeeze(max(CS_data.num_isolating_combined))/cmax;

    % Check against the threshold criteria and store
    heatmap_vals = sum(peak_isol_CS_data>0.05)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
end

if sum(strcmp(variable_name,'duration') == 1)

    % Peak isolations
    prev_data = CS_data.prevpresymp_combined + CS_data.prevasymp_combined + CS_data.prevsymp_combined;
    duration_data = zeros(size(squeeze(prev_data(1,:,:))));
    for i = 1:length(CS_data.numinf_combined(1,:,1))
        for j = 1:length(CS_data.numinf_combined(1,1,:))
            duration_data(i,j) = find(prev_data(:,i,j)>0,1,'last');        
        end
    end
    
    % Check against the threshold criteria and store
    heatmap_vals = sum(duration_data>150)/n_simns;

    % Reorder into a 2D array. Row by transrisk, column by worker group size
    heatmap_array = reshape(heatmap_vals,4,3);
end

end
