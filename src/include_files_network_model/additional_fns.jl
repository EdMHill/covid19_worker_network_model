"""
Purpose:
Houses functions required to run worker pattern network model

The functions contained within this file include:
- populate_atwork                   (arrays specifying when individuals are at the workplace)
- increment_counters!               (Used each timestep to increase the value of time in state variables)
- load_configs                      (For running batches of scenarios, selected by string value)
- find_network_parameters           (Load relevant network params based on number of worker types request)

Transmission functions contained within this file include:
- transmit_over!                    (infection event check over list of contacts, though not used for other workplace contacts)
- transmit_over_other_workplace!    (infection event check over list of other workplace contacts)

Functions to set up transmission rates within household for each individual
- assign_household_transmit_onegroup!  (use if everyone has the same household transmission risk)

Functions to reinitialise states at start of each run
- reinitialise_node_states!         (multiply time series vectors by 0)
- reinitialise_daily_record_arrays! (reset arrays logging isolation and work attendance by day)
- reinitialise_workplace_params!    (reinitialise the workplace related quantities)
- reinitialise_CT_vars!             (reinitialise the contact tracing variables)

Miscellaneous functions contained within this file include:
- draw_sample_from_pmf!
- set_infection_related_times!      (set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence)
"""

function populate_atwork!(atwork::Array{Int64,2},
                            sameday::Int64,
                            ton::Int64,
                            toff::Int64,
                            rng::MersenneTwister)


# Initialise array to store indicator values of whether individuals are at workplace or not each day
lmul!(0,atwork)

cmax,endtime = size(atwork)

# If sameday=0, all workers are at work on the same set of consecutive days.
    # If sameday=1, workers go to work on a random set of consecutive days.
    # If sameday=2, workers go to work on the same number of days, but scattered randomly throughout the week.
    if sameday==0 # all workers go to work on the same set of consecutive days
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                for days=ton:toff # all works work days ton to toff
                    atwork[node_itr,(reps_itr-1)*7+(days+1)] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            for days=ton:toff # all works work days ton to toff
                if (num_reps*7+(days+1))<=endtime
                    atwork[node_itr,num_reps*7+(days+1)] = 1
                end
            end
        end
    elseif sameday==1 # workers go to work on a random set of consecutive days

        # initialise pap
        pap = zeros(7)
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            dayon = ceil(Int64,7*rand(rng)) # which day do they start work
            lmul!(0,pap) # reinitalise pap
            # stay at work for toff days
            for t_it = 0:toff

                # Get day number of week.
                # uses mod1, rather than mod, so mod1(7,7) returns 7 rather than 0
                day_of_week = mod1(dayon+t_it,7)

                # Give value to pap
                pap[day_of_week] = 1
            end

            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the days at work in pap
                for day_itr = 1:7
                    atwork[node_itr,(reps_itr-1)*7+day_itr] = pap[day_itr]
                end
            end

            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for day_itr = 1:7
                if (num_reps*7+day_itr)<=endtime
                    atwork[node_itr,num_reps*7+day_itr] = pap[day_itr]
                end
            end
        end
    elseif sameday==2 # workdays randomly placed throughout the week (but repeat each week)
        # iterate over each node
        pap = zeros(Int64,7)
        for node_itr=1:cmax
            randperm!(rng,pap)  # Get permutation of 1:7
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the first toff days in the random permuation, pap
                for pap_itr = 1:(toff+1)
                    atwork[node_itr,(reps_itr-1)*7+pap[pap_itr]] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for pap_itr = 1:(toff+1)
                if (num_reps*7+pap[pap_itr])<=endtime
                    atwork[node_itr,num_reps*7+pap[pap_itr]] = 1
                end
            end
        end
    elseif sameday==3 # do ton weeks on followed by toff weeks off, all simultaneous
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if mod1(reps_itr,ton+toff) <= ton # using mod1, as mod1(1,1) = 1, whereas mod(1,1) = 0
                    for days=1:5 # only work weekdays
                        atwork[node_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    elseif sameday==4 # do ton weeks on followed by toff weeks off, starting at a random week each
        # initialise pap
        pap = zeros(ton+toff)
        # iterate over each node
        for node_itr=1:cmax
            # which week will they start
            weekon = ceil(Int64,(ton+toff)*rand(rng))
            lmul!(0,pap) # reinitalise pap
            # stay at work for ton weeks
            for t_it = 1:ton
                pap[mod1(weekon+(t_it-1),ton+toff)] = 1
            end
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if pap[mod1(reps_itr,ton+toff)] == 1
                    for days=1:5 # only work weekdays
                        atwork[node_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    end

    return nothing
end

function increment_counters!(states::node_states,
    household_isoltime::Int64, symp_isoltime::Int64,contact_tracing_active::Bool;
    timeisol_CTcause::Array{Int64,1}=zeros(Int64,1),CT_caused_isol_limit::Int64=0)

@unpack timelat, timeinf, timesymp, timeisol, symp_timeisol, lattime, inftime,
symptime, asymp, timeisol, symp_timeisol = states

# Increments time counters
    cmax = length(timelat)
    for node_itr = 1:cmax
        if timelat[node_itr]>0
            timelat[node_itr] += 1
        end

        if timeinf[node_itr]>0
            timeinf[node_itr] += 1
        end

        if timesymp[node_itr]>0
            timesymp[node_itr] += 1
        end

        if timeisol[node_itr]>0
            timeisol[node_itr] += 1
        end

        if timeisol[node_itr]>household_isoltime
            timeisol[node_itr] = 0
        end

        if symp_timeisol[node_itr]>0
            symp_timeisol[node_itr] += 1
        end

        if symp_timeisol[node_itr]>symp_isoltime
            symp_timeisol[node_itr] = 0
            timeisol[node_itr] = 0
            if contact_tracing_active == true
                timeisol_CTcause[node_itr] = 0
            end
        end

        if contact_tracing_active == true

            # Increment time in self-isolation, caused by Contact Tracing
            if timeisol_CTcause[node_itr]>0
                timeisol_CTcause[node_itr] += 1
            end

            # Reset contact tracing counter if limit is exceeded
            if timeisol_CTcause[node_itr]>CT_caused_isol_limit
                timeisol_CTcause[node_itr] = 0
            end
        end
    end
    return nothing
end


function load_configs(runset::String,workertypes::Int64,cmax::Int64,
                        RNGseed::Int64
                        )

# first set default parameters (then overwrite those that need to be changed)
toff = 4
ton = 0
sameday = 2
work_percent = 1.0

# workplaces close when 50% report infection
# only used if workplace_closure_active = true
workplace_CT_threshold = 0.5

# 10% of people correctly identify infector
# only used if perform_CT_from_infector = true
prob_backwards_CT = 0.
infector_engage_with_CT_prob = 1.0

# Scale isolated transmission risk
iso_trans_scaling = 1.

# (CT_engagement)% of those adhering to self isolation engage with contact tracing
# only used if contact_tracing_active = true
CT_engagement = 1.

# Set defaults for sensitivity configs
trans_scaling = [0.8,0.8,0.8,0.8,0.8]
adherence = 0.
clustering = [0.05, 0.5]
CS_team_size = 2

# Initialise empty intervention list
intervention_list = intervention_struct[]


# If sameday=0, all workers are at work on the same set of consecutive days.
# If sameday=1, workers go to work on a random set of consecutive days.
# If sameday=2, workers go to work on the same number of days, but scattered randomly throughout the week.
# If sameday=3, workers go to work ton weeks on followed by toff weeks off, all in synchrony
# If sameday=4, workers go to work ton weeks on followed by toff weeks off, beginning at a random time
if runset=="synchronised_changedays"
    toff_config = [-1:4;]
    n_configs = length(toff_config)
    sameday = 0
elseif runset=="variable_changedays"
    toff_config = [-1:4;]
    n_configs = length(toff_config)
elseif runset == "groups_off"
    if workertypes!=4
        display("error - for groups off, need 4 workertypes")
    else
        toff_config = ones(Int64,4)*4
        ton_config = zeros(Int64,length(toff_config))
        sameday_config = ones(Int64,length(toff_config))*0
        work_percent_config = [    [0  1  1  1]
            [1  0  1  1]
            [1  1  0  1]
            [1  1  1  0]]*0.5
    end
elseif runset=="weeks_on"
    toff_config = [1 1 2 2 2 2]
    ton_config = [1 1 2 2 1 1]
    sameday_config = [3 4 3 4 3 4]
    n_configs = length(toff_config)
elseif runset=="amount_backwards_CT"
    prob_backwards_CT_config = [0:0.05:0.5;]
    n_configs = length(prob_backwards_CT_config)
elseif runset=="run_one_run"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    adherence = 0.
    n_configs = 1
    transiso_config = ones(Float64, n_configs)
elseif runset=="workplace_CT_threshold"
    workplace_CT_threshold_config = [0:0.05:0.5;]
    n_configs = length(workplace_CT_threshold_config)
elseif runset=="CT_engagement"
    CT_engagement_config = [0:0.1:1;]
    n_configs = length(CT_engagement_config)
elseif runset=="dynamic_social_timeframe"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 0.5
    dynamic_time_frame_config = [1,2,3,5,7,10,14]
    n_configs = length(dynamic_time_frame_config)
elseif runset=="dynamic_social_timeframe_groups"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 0.5
    dynamic_time_frame_config = [1,2,3,5,7,10,14]
    n_configs = length(dynamic_time_frame_config)
elseif runset=="social_group_size"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    max_contacts_social_config = [3,6,12,24,48,100]
    n_configs = length(max_contacts_social_config)
elseif runset=="rule_of_n"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    group_limit_config = [2,3,4,6,12]
    n_configs = length(group_limit_config)
elseif runset=="rule_of_n_weekly"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    group_limit_config = [2,3,4,6,12]
    n_configs = length(group_limit_config)
elseif runset=="rule_of_n_monthly"
    # everyone works mon-fri
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    group_limit_config = [2,3,4,6,12]
    n_configs = length(group_limit_config)
elseif runset=="lockdown_adherence"
elseif runset=="change_cmax"
    cmax_config = [10000:10000:50000;]
    n_configs = length(cmax_config)
elseif runset=="test_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    intervention1 = intervention_struct(effects = ["workpercent"],
                                        start_time=30,
                                        workpercent=zeros(workertypes))
    push!(intervention_list, intervention1)
    intervention2 = intervention_struct(effects = ["workpercent"],
                                        start_time = 130,
                                        workpercent=ones(workertypes))
    push!(intervention_list, intervention2)
elseif runset=="lockdown_duration_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    duration_intervention_options = [20:20:80;]
    n_int_sets = length(duration_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["workpercent"],
                                        start_time=30,
                                        workpercent=zeros(workertypes)),
                            intervention_struct(effects = ["workpercent"],
                                        start_time = 30+duration_intervention_options[i],
                                        workpercent=ones(workertypes))] for i = 1:n_int_sets]

elseif runset=="synchronised_changedays_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    toff_intervention_options = [-1:4;]
    n_int_sets = length(toff_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["worker_patterns","adherence","contact_tracing"],
                                            start_time = 15,
                                            sameday = 0,
                                            toff = toff_intervention_options[i],
                                            ton = 0,
                                            CT_parameters = CT_params(contact_tracing_active=true,
                                                                    CT_engagement=1.))] for i = 1:n_int_sets]
elseif runset=="variable_changedays_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    toff_intervention_options = [-1:4;]
    n_int_sets = length(toff_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["worker_patterns","adherence","contact_tracing"],
                                            start_time = 15,
                                            sameday = 2,
                                            toff = toff_intervention_options[i],
                                            ton = 0,
                                            CT_parameters = CT_params(contact_tracing_active=true,
                                                                    CT_engagement=1.))] for i = 1:n_int_sets]
elseif runset=="weeks_on_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    toff_intervention_options = [1 1 2 2 2 2]
    ton_intervention_options = [1 1 2 2 1 1]
    sameday_intervention_options = [3 4 3 4 3 4]
    n_int_sets = length(toff_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["worker_patterns"],
                                            start_time = 30,
                                            sameday = sameday_intervention_options[i],
                                            toff = toff_intervention_options[i],
                                            ton = ton_intervention_options[i])] for i = 1:n_int_sets]
elseif runset=="workpercent_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workpercent_intervention_options = [1:-0.1:0;]
    n_int_sets = length(workpercent_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["workpercent","adherence","contact_tracing"],
                                            start_time = 15,
                                            workpercent = workpercent_intervention_options[i]*ones(workertypes),
                                            CT_parameters = CT_params(contact_tracing_active=true,
                                                                    CT_engagement=1.))
                                            ] for i = 1:n_int_sets]

    # Add scenario where proportion of each type of worker returning to work is not constant across all sectors
    # work_percent_config_segment2 = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]
    work_percent_config_segment2 = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]
    intervention_config2 = intervention_struct(effects = ["workpercent","adherence","contact_tracing"],
                                            start_time = 15,
                                            workpercent = work_percent_config_segment2,
                                            CT_parameters = CT_params(contact_tracing_active=true,
                                                                    CT_engagement=1.))
    # Concatenate the two sets of configs
    push!(intervention_list_config,[intervention_config2])
elseif runset=="amount_backwards_CT_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    prob_backwards_CT_intervention_options = [0:0.1:1.0;]
    n_int_sets = length(prob_backwards_CT_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["contact_tracing","adherence"],
                                                    start_time = 15,
                                                    CT_parameters = CT_params(contact_tracing_active=true,
                                                                        prob_backwards_CT=prob_backwards_CT_intervention_options[i],
                                                                        CT_engagement=CT_engagement,
                                                                        perform_CT_from_infector=true))]
                                                                        for i = 1:n_int_sets]
elseif runset=="adherence_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    adherence_intervention_options = [0:0.1:1;]
    n_int_sets = length(adherence_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["adherence","contact_tracing"],
                                            start_time = 15,
                                            adherence = adherence_intervention_options[i],
                                            CT_parameters = CT_params(contact_tracing_active=true,
                                                                    CT_engagement=1.))]
                                            for i = 1:n_int_sets]
elseif runset=="CS_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    CS_team_size_options = [2, 5, 10]
    CS_scale_transrisk_options = [0.25,0.5,0.75,1]
    n_CS_team_size_ops = length(CS_team_size_options)
    n_CS_scale_transrisk_ops = length(CS_scale_transrisk_options)
    n_int_sets = n_CS_team_size_ops*n_CS_scale_transrisk_ops
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)
    set_idx = 1

    for CS_team_size_itr = 1:n_CS_team_size_ops
        CS_team_size_interv = CS_team_size_options[CS_team_size_itr]
        for CS_scale_transrisk_itr = 1:n_CS_scale_transrisk_ops
            CS_scale = CS_scale_transrisk_options[CS_scale_transrisk_itr]

            # Set up intervention list
            intervention_list_config[set_idx] = [intervention_struct(effects = ["COVID_secure","adherence","contact_tracing"],
                                                    start_time = 15,
                                                    CS_scale_transrisk_scalar = CS_scale,
                                                    CS_team_size_intervention = CS_team_size_interv,
                                                    CT_parameters = CT_params(contact_tracing_active=true,
                                                                            CT_engagement=1.))]

            # Increment set index
            set_idx += 1
        end
    end
elseif runset=="CS_intervention_no_isol"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    CS_team_size_options = [2, 5, 10]
    CS_scale_transrisk_options = [0.25,0.5,0.75,1]
    n_CS_team_size_ops = length(CS_team_size_options)
    n_CS_scale_transrisk_ops = length(CS_scale_transrisk_options)
    n_int_sets = n_CS_team_size_ops*n_CS_scale_transrisk_ops
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)
    set_idx = 1

    for CS_team_size_itr = 1:n_CS_team_size_ops
        CS_team_size_interv = CS_team_size_options[CS_team_size_itr]
        for CS_scale_transrisk_itr = 1:n_CS_scale_transrisk_ops
            CS_scale = CS_scale_transrisk_options[CS_scale_transrisk_itr]

            # Set up intervention list
            intervention_list_config[set_idx] = [intervention_struct(effects = ["COVID_secure"],
                                                    start_time = 15,
                                                    CS_scale_transrisk_scalar = CS_scale,
                                                    CS_team_size_intervention = CS_team_size_interv)]

            # Increment set index
            set_idx += 1
        end
    end
elseif runset=="CS_workplace_no_control"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [2, 5, 10, 15, 25, 50, 100, 1000]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    CS_team_size_config = variable_ops
elseif runset=="workplace_CT_threshold_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workplace_CT_threshold_intervention_options = [0:0.05:0.4;]
    n_int_sets = length(workplace_CT_threshold_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["contact_tracing"],
                                                    start_time = 30,
                                                    CT_parameters = CT_params(contact_tracing_active=true,
                                                                        workplace_CT_threshold=workplace_CT_threshold_intervention_options[i],
                                                                        CT_engagement=CT_engagement,
                                                                        workplace_closure_active=true))]
                                                                        for i = 1:n_int_sets]
elseif runset=="hardfast_vs_weakslow_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workpercent_intervention_options = [0.2:0.2:0.8;]
    intervention_length_options = [40:40:160;]
    n_int_sets = length(workpercent_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["workpercent"],
                                            start_time = 30,
                                            workpercent = workpercent_intervention_options[i]*ones(workertypes)),
                                intervention_struct(effects = ["workpercent"],
                                                    start_time = 30+intervention_length_options[i],
                                                    workpercent = ones(workertypes))] for i = 1:n_int_sets]
elseif runset=="adaptive_workplace_closures"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    max_lockdown_weeks = 4
    workpercent_intervention_options = [0, 0.5, 0.75]
    n_blocks_intervention_options = [1, 2, 4]
    block_gap_weeks_options = [1, 2, 4]          # weeks between lockdowns
    n_int_sets = length(workpercent_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)
    workpercent_op_count = 1
    n_block_count = 1
    block_gap_count = 1
    for int_set_itr = 1:n_int_sets
        workpercent = workpercent_intervention_options[workpercent_op_count]
        n_blocks = n_blocks_intervention_options[n_block_count]
        block_gap = block_gap_weeks_options[block_gap_count]*7

        block_length = Int64((max_lockdown_weeks/(1-workpercent))/n_blocks)*7

        intervention_list = [intervention_struct(effects = ["workpercent"],
                                            start_time = 30,
                                            workpercent = workpercent*ones(workertypes)),
                            intervention_struct(effects = ["workpercent"],
                                            start_time = 30+block_length,
                                            workpercent = ones(workertypes))]

        for block_itr = 1:(n_blocks-1)
            push!(intervention_list, intervention_struct(effects = ["workpercent"],
                                                start_time = 30+block_itr*(block_length+block_gap),
                                                workpercent = workpercent*ones(workertypes)))
            push!(intervention_list, intervention_struct(effects = ["workpercent"],
                                                start_time = 30+block_itr*(block_length+block_gap)+block_length,
                                                workpercent = ones(workertypes)))
        end
        intervention_list_config[int_set_itr] = intervention_list

        block_gap_count += 1
        if block_gap_count > length(block_gap_weeks_options)
            block_gap_count = 1
            n_block_count += 1
        end
        if n_block_count > length(n_blocks_intervention_options)
            n_block_count = 1
            workpercent_op_count += 1
        end
    end

elseif runset=="multiphase_lockdown"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workpercent_intervention_options = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    start_intervention_options = [30, 60, 90, 120]
    length_intervention_options = [1, 2, 4, 8, 12, 16]
    n_blocks_intervention_options = [1, 2, 4, 6]
    block_gap_weeks_options = [1, 2, 4, 8, 12]          # weeks between lockdowns
    n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)
    workpercent_op_count = 1
    start_op_count = 1
    length_op_count = 1
    n_block_count = 1
    block_gap_count = 1
    for int_set_itr = 1:n_int_sets
        workpercent = workpercent_intervention_options[workpercent_op_count]
        start_time = start_intervention_options[start_op_count]
        total_length = length_intervention_options[length_op_count]
        n_blocks = n_blocks_intervention_options[n_block_count]
        block_gap = block_gap_weeks_options[block_gap_count]*7

        block_length = floor(Int64, (total_length*7)/n_blocks)

        intervention_list = [intervention_struct(effects = ["workpercent","transrisk"],
                                            start_time = start_time,
                                            workpercent = workpercent*ones(workertypes),
                                            scaling_social = workpercent,
                                            scaling_random = workpercent),
                            intervention_struct(effects = ["workpercent","transrisk"],
                                            start_time = start_time+block_length,
                                            workpercent = ones(workertypes))]

        for block_itr = 1:(n_blocks-1)
            push!(intervention_list, intervention_struct(effects = ["workpercent","transrisk"],
                                                start_time = start_time+block_itr*(block_length+block_gap),
                                                workpercent = workpercent*ones(workertypes),
                                                scaling_social = workpercent,
                                                scaling_random = workpercent))
            push!(intervention_list, intervention_struct(effects = ["workpercent","transrisk"],
                                                start_time = start_time+block_itr*(block_length+block_gap)+block_length,
                                                workpercent = ones(workertypes)))
        end
        intervention_list_config[int_set_itr] = intervention_list

        block_gap_count += 1
        if block_gap_count > length(block_gap_weeks_options)
            block_gap_count = 1
            n_block_count += 1
        end
        if n_block_count > length(n_blocks_intervention_options)
            n_block_count = 1
            length_op_count += 1
        end
        if length_op_count > length(length_intervention_options)
            length_op_count = 1
            start_op_count += 1
        end
        if start_op_count > length(start_intervention_options)
            start_op_count = 1
            workpercent_op_count += 1
        end
    end

elseif runset=="multiphase_lockdown_contact_struct"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workpercent_intervention_options = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    start_intervention_options = [30, 60, 90, 120]
    length_intervention_options = [1, 2, 4, 8, 12, 16]
    n_blocks_intervention_options = [1, 2, 4, 6]
    block_gap_weeks_options = [1, 2, 4, 8, 12]          # weeks between lockdowns
    n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)
    workpercent_op_count = 1
    start_op_count = 1
    length_op_count = 1
    n_block_count = 1
    block_gap_count = 1
    for int_set_itr = 1:n_int_sets
        workpercent = workpercent_intervention_options[workpercent_op_count]
        start_time = start_intervention_options[start_op_count]
        total_length = length_intervention_options[length_op_count]
        n_blocks = n_blocks_intervention_options[n_block_count]
        block_gap = block_gap_weeks_options[block_gap_count]*7

        block_length = floor(Int64, (total_length*7)/n_blocks)

        intervention_list = [intervention_struct(effects = ["workpercent","contact_structure"],
                                            start_time = start_time,
                                            workpercent = workpercent*ones(workertypes),
                                            social_contact_scaling = workpercent,
                                            random_contact_scaling = workpercent),
                            intervention_struct(effects = ["workpercent","contact_structure"],
                                            start_time = start_time+block_length,
                                            workpercent = ones(workertypes))]

        for block_itr = 1:(n_blocks-1)
            push!(intervention_list, intervention_struct(effects = ["workpercent","contact_structure"],
                                                start_time = start_time+block_itr*(block_length+block_gap),
                                                workpercent = workpercent*ones(workertypes),
                                                social_contact_scaling = workpercent,
                                                random_contact_scaling = workpercent))
            push!(intervention_list, intervention_struct(effects = ["workpercent","contact_structure"],
                                                start_time = start_time+block_itr*(block_length+block_gap)+block_length,
                                                workpercent = ones(workertypes)))
        end
        intervention_list_config[int_set_itr] = intervention_list

        block_gap_count += 1
        if block_gap_count > length(block_gap_weeks_options)
            block_gap_count = 1
            n_block_count += 1
        end
        if n_block_count > length(n_blocks_intervention_options)
            n_block_count = 1
            length_op_count += 1
        end
        if length_op_count > length(length_intervention_options)
            length_op_count = 1
            start_op_count += 1
        end
        if start_op_count > length(start_intervention_options)
            start_op_count = 1
            workpercent_op_count += 1
        end
    end
elseif runset=="transrisk_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    scaling_intervention_options = [0:0.2:1;]
    n_int_sets = length(scaling_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["transrisk"],
                                            start_time = 30,
                                            scaling_work_static = scaling_intervention_options[i],
                                            scaling_work_dynamic = scaling_intervention_options[i],
                                            scaling_social = scaling_intervention_options[i],
                                            scaling_random = scaling_intervention_options[i]),
                                intervention_struct(effects = ["transrisk"],
                                            start_time = 80)
                                            ] for i = 1:n_int_sets]

elseif runset=="CT_engagement_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    CT_engagement_intervention_options = [0:0.2:1;]
    n_int_sets = length(CT_engagement_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["contact_tracing"],
                                                    start_time = 30,
                                                    CT_parameters = CT_params(contact_tracing_active=true,
                                                                            CT_engagement=CT_engagement_intervention_options[i]))]
                                                                            for i = 1:n_int_sets]
elseif runset=="contact_structure_intervention"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    contact_scaling_intervention_options = [0:0.2:1;]
    n_int_sets = length(contact_scaling_intervention_options)
    intervention_list_config = [[intervention_struct(effects = ["contact_structure"],
                                            start_time = 30,
                                            social_contact_scaling = contact_scaling_intervention_options[i],
                                            random_contact_scaling = contact_scaling_intervention_options[i]),
                                intervention_struct(effects = ["contact_structure"],
                                            start_time = 80)
                                            ] for i = 1:n_int_sets]

elseif runset=="adaptive_lockdown_adherence"
    n_configs = 1
    sameday = 3
    ton = 1
    toff = 0
    workpercent_intervention_options = [0, 0.5, 1.0]
    start_intervention_options = [30, 60, 90]
    length_intervention_options = [1, 2, 4, 8]
    n_blocks_intervention_options = [1, 2]
    block_gap_weeks_options = [1, 2, 4]          # weeks between lockdowns
    n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
    intervention_list_config = Array{Array{intervention_struct,1},1}(undef, n_int_sets)

    adherence_period = 28   # time scale (days) at which adherence naturally decreases
    n_adherence_periods = ceil(Int64, endtime/adherence_period)
    adherence_period_decrease = 0.05

    workpercent_op_count = 1
    start_op_count = 1
    length_op_count = 1
    n_block_count = 1
    block_gap_count = 1
    for int_set_itr = 1:n_int_sets
        workpercent = workpercent_intervention_options[workpercent_op_count]
        start_time = start_intervention_options[start_op_count]
        total_length = length_intervention_options[length_op_count]
        n_blocks = n_blocks_intervention_options[n_block_count]
        block_gap = block_gap_weeks_options[block_gap_count]*7

        block_length = floor(Int64, (total_length*7)/n_blocks)

        policy_change_times = zeros(Int64, n_blocks*2)
        lockdown_active_array = zeros(Int64, endtime)
        for block_itr = 1:n_blocks
            policy_change_times[(2*(block_itr-1)+1)] = start_time+(block_itr-1)*(block_length+block_gap)
            policy_change_times[(2*block_itr)] = start_time+(block_itr-1)*(block_length+block_gap)+block_length
            lockdown_active_array[(start_time+(block_itr-1)*(block_length+block_gap)):(start_time+(block_itr-1)*(block_length+block_gap)+block_length)] .= 1
        end

        cumulative_lockdown_days = cumsum(lockdown_active_array)
        cumulative_lockdown_periods = floor.(Int64, cumulative_lockdown_days./adherence_period)

        adherence_over_time = ones(Float64, endtime).*adherence
        for time_itr = 1:endtime
            adherence_over_time[time_itr] = adherence*(1-adherence_period_decrease)^cumulative_lockdown_periods[time_itr]
        end

        intervention_list = intervention_struct[]

        time_itr = adherence_period + 1
        while time_itr <= endtime

            if (time_itr+adherence_period-1) <= endtime
                current_time_period = [time_itr:(time_itr+adherence_period-1);]
            else
                current_time_period = [time_itr:endtime;]
            end

            # Find policy changes in current time period
            policy_change_ids = findall((policy_change_times.>=current_time_period[1]).*(policy_change_times.<=current_time_period[end]).==1)
            n_policy_change = length(policy_change_ids)

            if adherence_over_time[time_itr] != adherence_over_time[time_itr-adherence_period]
                push!(intervention_list, intervention_struct(effects = ["adherence"],
                                                start_time = time_itr,
                                                adherence = adherence_over_time[time_itr]))
            end
            if n_policy_change > 0
                for policy_change_itr = 1:n_policy_change
                    push!(intervention_list, intervention_struct(effects = ["workpercent","contact_structure"],
                                                    start_time = policy_change_times[policy_change_ids[policy_change_itr]],
                                                    workpercent = workpercent*ones(workertypes),
                                                    social_contact_scaling = workpercent,
                                                    random_contact_scaling = workpercent))
                end

            end
            time_itr += adherence_period
        end

        intervention_list_config[int_set_itr] = intervention_list

        block_gap_count += 1
        if block_gap_count > length(block_gap_weeks_options)
            block_gap_count = 1
            n_block_count += 1
        end
        if n_block_count > length(n_blocks_intervention_options)
            n_block_count = 1
            length_op_count += 1
        end
        if length_op_count > length(length_intervention_options)
            length_op_count = 1
            start_op_count += 1
        end
        if start_op_count > length(start_intervention_options)
            start_op_count = 1
            workpercent_op_count += 1
        end
    end

#### Sensitivity configs, assuming no intervention or isolation ####
elseif runset=="RNGseed_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [100, 200, 300]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    rng_config = variable_ops
    adherence = 0.
# elseif runset=="initinfected_svty"
#     sameday = 3
#     ton = 1
#     toff = 0
#     work_percent = 1.0
#     variable_ops = [1, 5, 10, 50]
#     n_configs = length(variable_ops)
#     transiso_config = ones(Float64, n_configs)
#     init_infected_config = variable_ops
elseif runset=="transscaling_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    adherence = 0.
    variable_ops = [[0.6 0.6 0.6 0.6 0.6], [0.7 0.7 0.7 0.7 0.7], [0.8 0.8 0.8 0.8 0.8],
                    [0.9 0.9 0.9 0.9 0.9], [1 1 1 1 1]]
    # variable_ops = [[0 1 1 1], [0.25 1 1 1], [0.5 1 1 1], [0.75 1 1 1],
    #                 [1 0 1 1], [1 0.25 1 1], [1 0.5 1 1], [1 0.75 1 1],
    #                 [1 1 0 1], [1 1 0.25 1], [1 1 0.5 1], [1 1 0.75 1],
    #                 [1 1 1 0], [1 1 1 0.25], [1 1 1 0.5], [1 1 1 0.75],
    #                 [0.3 0.3 0.3 0.3], [0.4 0.4 0.4 0.4], [0.5 0.5 0.5 0.5], [0.6 0.6 0.6 0.6], [0.7 0.7 0.7 0.7], [0.8 0.8 0.8 0.8], [0.9 0.9 0.9 0.9], [1 1 1 1 1], [1.1 1.1 1.1 1.1], [1.2 1.2 1.2 1.2], [1.3 1.3 1.3 1.3], [1.4 1.4 1.4 1.4], [1.5 1.5 1.5 1.5], [1.6 1.6 1.6 1.6], [1.7 1.7 1.7 1.7], [1.8 1.8 1.8 1.8], [1.9 1.9 1.9 1.9], [2 2 2 2]]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    trans_scaling_config = variable_ops
elseif runset=="adherence_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    adherence_config = variable_ops
elseif runset=="clustering_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [[0.05, 0.5]]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="social_clustering_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.25, 0.5, 0.75, 1]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="popsize_svty"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [5000, 10000, 25000, 50000, 100000]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    cmax_config = variable_ops
elseif runset=="ER_no_control"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    n_configs = 1
    transiso_config = ones(Float64, n_configs)
elseif runset=="CS_workplace_no_control"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [2, 5, 10, 15, 25, 50, 100, 1000]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    CS_team_size_config = variable_ops
elseif runset=="workplace_clustering_svty_no_social"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_no_household"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_workplace_only"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_no_dynamic"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_workplace_static_only"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_workplace_static_and_social"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="workplace_clustering_svty_workplace_static_and_household"
    sameday = 3
    ton = 1
    toff = 0
    work_percent = 1.0
    variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
    n_configs = length(variable_ops)
    transiso_config = ones(Float64, n_configs)
    clustering_config = variable_ops
elseif runset=="work_percent_svty"
    sameday = 3
    ton = 1
    toff = 0
    variable_ops = [0.2,0.4,0.6,0.8,1.0]
    n_configs = length(variable_ops) + 1
    transiso_config = ones(Float64, n_configs)
    work_percent_config_segment1 = repeat(variable_ops, outer=[1,workertypes])
    work_percent_config_segment2 = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]   # proportion of each type of worker returning to work
    work_percent_config = vcat(work_percent_config_segment1,work_percent_config_segment2')
end
##### End of sensitivity configs #####

# extend scalars to vectors where needed
if @isdefined(sameday_config)==false
    sameday_config = sameday*ones(Int64,n_configs)
end
if @isdefined(ton_config)==false
    ton_config = ton*ones(Int64,n_configs)
end
if @isdefined(toff_config)==false
    toff_config = toff*ones(Int64,n_configs)
end
if @isdefined(work_percent_config)==false
    work_percent_config = work_percent*ones(n_configs,workertypes)
end
if @isdefined(workplace_CT_threshold_config)==false
    workplace_CT_threshold_config = workplace_CT_threshold*ones(n_configs)
end
if @isdefined(prob_backwards_CT_config)==false
    prob_backwards_CT_config = prob_backwards_CT*ones(n_configs)
end
if @isdefined(infector_engage_with_CT_prob_config)==false
    infector_engage_with_CT_prob_config = infector_engage_with_CT_prob*ones(n_configs)
end
if @isdefined(CT_engagement_config)==false
    CT_engagement_config = CT_engagement*ones(n_configs)
end
# if @isdefined(transasymp_config)==false
#     transasymp_config = asymp_trans_scaling*ones(n_configs)
# end
if @isdefined(transiso_config)==false
    transiso_config = iso_trans_scaling*ones(n_configs)
end
if @isdefined(cmax_config)==false
    cmax_config = cmax*ones(Int64,n_configs)
end
if @isdefined(rng_config)==false
    rng_config = RNGseed*ones(Int64,n_configs)
end
if @isdefined(trans_scaling_config)==false
    trans_scaling_config = [trans_scaling for i=1:n_configs]
end
# if @isdefined(init_infected_config)==false
#     # Specify number of symptomatic infecteds
#     init_infected_config = n_initial_infecteds[3]*ones(Int64,n_configs)
# end
if @isdefined(adherence_config)==false
    adherence_config = adherence*ones(Int64,n_configs)
end
if @isdefined(clustering_config)==false
    clustering_config = [clustering for i=1:n_configs]
end
if @isdefined(CS_team_size_config)==false
    CS_team_size_config = [CS_team_size for i=1:n_configs]
end
if @isdefined(intervention_list_config)==false
    intervention_list_config = [intervention_list]
end
if @isdefined(dynamic_time_frame_config)==false
    dynamic_time_frame_config = []
end
if @isdefined(group_limit_config)==false
    group_limit_config = []
end
if @isdefined(max_contacts_social_config)==false
    max_contacts_social_config = []
end

#### NEW CONFIGS ADDED HERE

    return sameday_config, ton_config, toff_config, work_percent_config,
    n_configs, workplace_CT_threshold_config, prob_backwards_CT_config,
    infector_engage_with_CT_prob_config, CT_engagement_config,
    transiso_config, cmax_config,
    rng_config, trans_scaling_config, adherence_config,
    clustering_config, CS_team_size_config, intervention_list_config,
    dynamic_time_frame_config, group_limit_config, max_contacts_social_config   #### NEW CONFIGS ADDED HERE
end

function find_network_parameters(workertypes;workpercent::Array{Float64,1}=Array{Float64,1}[])

    if workertypes==2
        if isempty(workpercent)
            workpercent = [0.5, 0.5]  # proportion of each type of worker returning to work
        end
        network_parameters = network_params(n_nodes = cmax,
            prob_workertype_contact = [0.002,0.002].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [10., 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [2., 0])            # sd for number of dynamic contacts for each worker type
        workplace_generation_parameters = workplace_generation_params(workertypes=2,
            workpercent=workpercent,
            workforce_proportion = [0.15, 0.85],   # proportion of workforce in each worker type
            workplace_size_mean = [10., 100.],      # mean size of workplace for each worker group
            workplace_size_sd = [2., 10. ])       # sd for size of workplace for each worker group
    elseif workertypes==4
        if isempty(workpercent)
            workpercent = [1.,1.,1.,1.]              # proportion of each type of worker returning to work
        end
        network_parameters = network_params(n_nodes = cmax,
            prob_workertype_contact = [0.0001,0.0001,0.0001,0.0001].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [10., 10., 10., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0., 0., 0., 10.],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0., 0., 0., 2.])            # sd for number of dynamic contacts for each worker type
        workplace_generation_parameters = workplace_generation_params(workertypes=4,
            workpercent=workpercent,
            workforce_proportion = [0.030, 0.127, 0.311, 0.532],   # proportion of workforce in each worker type
            workplace_size_mean = [6., 9., 11., 10.],      # mean size of workplace for each worker group
            workplace_size_sd = [14., 57., 122., 109. ])       # sd for size of workplace for each worker group
    elseif workertypes==6
        if isempty(workpercent)
            workpercent = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]            # proportion of each type of worker returning to work
        end
        network_parameters = network_params(n_nodes = cmax,
            prob_workertype_contact = [0.002,0.002,0.002,0.002,0.002,0.002].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 20., 20.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0., 0., 0., 10., 10., 10.],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0., 0., 0., 2., 2., 2.])            # sd for number of dynamic contacts for each worker type
        workplace_generation_parameters = workplace_generation_params(workertypes=6,
            workpercent=workpercent,
            workforce_proportion = [0.025, 0.105, 0.252, 0.427, 0.072, 0.119],   # proportion of workforce in each worker type
            workplace_size_mean = [6., 9., 11., 10., 65., 33.],      # mean size of workplace for each worker group
            workplace_size_sd = [14., 57., 122., 109., 432., 305. ])       # sd for size of workplace for each worker group
    elseif workertypes==41
        if isempty(workpercent)
            workpercent = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]            # proportion of each type of worker returning to work
        end
        network_parameters = network_params(n_nodes = cmax,
            prob_workertype_contact = ones(workertypes)*0.002.*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 10., 10., 10., 10., 10., 10., 10., 10., 5., 5., 10., 5., 5., 5., 10., 10., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 5., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0, 0, 0, 0, 2., 0, 5., 2., 10., 10., 0, 10., 10., 10., 0, 0, 0, 0, 0, 0, 5., 5., 0, 5., 0, 2., 0, 0, 2., 10., 2., 5., 10., 5., 10., 10., 0, 5., 5., 5., 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0, 0, 0, 0, 1., 0, 2., 1., 4., 4., 0, 4., 4., 4., 0, 0, 0, 0, 0, 0, 2., 2., 0, 2., 0, 1., 0, 0, 1., 4., 1., 2., 4., 2., 4., 4., 0, 2., 2., 2., 0],            # sd for number of dynamic contacts for each worker type
            workplace_degree_distribution = find_fitted_contact_dds("workplace_fixed",workertypes),
            workplace_dynamic_degree_distribution = find_fitted_contact_dds("workplace_dynamic",workertypes))
            workplace_generation_parameters = workplace_generation_params(workertypes=workertypes,
            workpercent=workpercent,
            workforce_proportion = [0.0393, 0.0035, 0.0181, 0.0958, 0.0126, 0.0627, 0.0263, 0.0582, 0.0700, 0.0206, 0.0124, 0.0032, 0.0218, 0.0674, 0.0137, 0.0273, 0.0028, 0.0188, 0.0207, 0.0957, 0.0021, 0.0063, 0.0260, 0.0033, 0.0044, 0.0184, 0.0194, 0.0064, 0.0614, 0.0399, 0.0306, 0.0342, 0.0078, 0.0016, 0.0144, 0.0027, 0.0131, 0.0039, 0.0074, 0.0010, 0.0049],   # proportion of workforce in each worker type
            workplace_size_mean = [5.3906, 57.6168, 34.7701, 15.5182, 19.2749, 3.7413, 6.9799, 11.5829, 6.8789, 6.1270, 12.8666, 2.9083, 24.7449, 9.9487, 5.6432, 3.3376, 6.4965, 6.3643, 4.2314, 4.1993, 10.8347, 7.2025, 16.9858, 8.0033, 10.1030, 8.6935, 3.3448, 17.3884, 28.3144, 14.6889, 59.8724, 19.5898, 5.0469, 31.3728, 10.4981, 8.1710, 12.1034, 5.7949, 3.5018, 10.8615, 4.0502],      # mean size of workplace for each worker group
            workplace_size_sd = [35.6780, 291.5813, 183.5026, 99.9800, 129.5788, 51.6139, 72.9294, 94.4625, 68.6747, 55.8071, 96.7104, 44.4382, 155.5328, 84.0205, 62.4778, 37.7083, 67.0943, 44.9687, 40.6361, 49.4326, 62.3247, 60.7611, 98.3582, 57.2267, 55.1600, 65.2562, 33.1987, 57.6397, 97.3243, 108.6843, 210.1163, 105.7828, 52.6915, 141.2081, 71.0224, 45.4627, 123.3399, 72.5046, 27.0021, 80.4898, 41.5465])       # sd for size of workplace for each worker group
    elseif workertypes==42
        if isempty(workpercent)
            workpercent = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3, 0.3]            # proportion of each type of worker returning to work
        end
        network_parameters = network_params(n_nodes = cmax,
            prob_workertype_contact = ones(workertypes)*0.002.*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 10., 10., 10., 10., 10., 10., 10., 10., 5., 5., 10., 5., 5., 5., 10., 10., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 5., 10., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0, 0, 0, 0, 2., 0, 5., 2., 10., 10., 0, 10., 10., 10., 0, 0, 0, 0, 0, 0, 5., 5., 0, 5., 0, 2., 0, 0, 2., 10., 2., 5., 10., 5., 10., 10., 0, 5., 5., 5., 0, 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0, 0, 0, 0, 1., 0, 2., 1., 4., 4., 0, 4., 4., 4., 0, 0, 0, 0, 0, 0, 2., 2., 0, 2., 0, 1., 0, 0, 1., 4., 1., 2., 4., 2., 4., 4., 0, 2., 2., 2., 0, 0])            # sd for number of dynamic contacts for each worker type
        workplace_generation_parameters = workplace_generation_params(workertypes=workertypes,
            workpercent=workpercent,
            workforce_proportion = [0.0393, 0.0035, 0.0181, 0.0958, 0.0126, 0.0627, 0.0263, 0.0582, 0.0700, 0.0206, 0.0124, 0.0032, 0.0218, 0.0674, 0.0137, 0.0273, 0.0028, 0.0188, 0.0207, 0.0957, 0.0021, 0.0063, 0.0260, 0.0033, 0.0044, 0.0184, 0.0194, 0.0064, 0.0614, 0.0399, 0.0306, 0.0342, 0.0078, 0.0016, 0.0144, 0.0027, 0.0131, 0.0039, 0.0074, 0.0010, 0.0049, 0.0001],   # proportion of workforce in each worker type
            workplace_size_mean = [5.3906, 57.6168, 34.7701, 15.5182, 19.2749, 3.7413, 6.9799, 11.5829, 6.8789, 6.1270, 12.8666, 2.9083, 24.7449, 9.9487, 5.6432, 3.3376, 6.4965, 6.3643, 4.2314, 4.1993, 10.8347, 7.2025, 16.9858, 8.0033, 10.1030, 8.6935, 3.3448, 17.3884, 28.3144, 14.6889, 59.8724, 19.5898, 5.0469, 31.3728, 10.4981, 8.1710, 12.1034, 5.7949, 3.5018, 10.8615, 4.0502, 79.6000],      # mean size of workplace for each worker group
            workplace_size_sd = [35.6780, 291.5813, 183.5026, 99.9800, 129.5788, 51.6139, 72.9294, 94.4625, 68.6747, 55.8071, 96.7104, 44.4382, 155.5328, 84.0205, 62.4778, 37.7083, 67.0943, 44.9687, 40.6361, 49.4326, 62.3247, 60.7611, 98.3582, 57.2267, 55.1600, 65.2562, 33.1987, 57.6397, 97.3243, 108.6843, 210.1163, 105.7828, 52.6915, 141.2081, 71.0224, 45.4627, 123.3399, 72.5046, 27.0021, 80.4898, 41.5465, 84.8072])       # sd for size of workplace for each worker group
    end

    return network_parameters, workplace_generation_parameters
end

function find_fitted_contact_dds(context::String, workertypes::Int64)

    if context == "workplace_fixed"
        if workertypes == 41
            work_contacts_meanlog = [1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.69, 1.69, 1.69, 1.67, 1.67, 3.48, 1.97, 2.86, 1.73, 1.73, 1.73, 3.1, 1.56, 1.5, 1.98, 1.56, 1.73, 1.97, 1.67, 1.35, 1.73, 1.67, 3.15, 1.96, 1.96, 1.96, 1.67, 1.67, 1.11, 1.67, 1.67, 1.73, 1.5, 1.97, 1.73]
            work_contacts_sdlog = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.2, 1.2, 1.2, 0.97, 0.97, 0.53, 1.35, 1.78, 0.91, 0.91, 0.91, 1.44, 0.33, 1.02, 0.92, 0.33, 0.91, 1.35, 1.04, 1.35, 0.91, 1.04, 1.43, 1.25, 1.25, 1.25, 0.49, 0.49, 1.12, 0.49, 1.04, 1.09, 0.11, 1.35, 0.91]

            workplace_degree_distribution = Array{Distribution,1}(undef, length(work_contacts_meanlog))

            for ii = 1:length(workplace_degree_distribution)
                workplace_degree_distribution[ii] = Distributions.LogNormal(work_contacts_meanlog[ii],work_contacts_sdlog[ii])
            end
        end

        return workplace_degree_distribution

    elseif context == "workplace_dynamic"
        if workertypes == 41
            work_contacts_meanlog = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 2.4, 2.4, 2.4, 3.05, 3.05, 1.76, 1.76, 1.76, 0.96, 0.96, 0.96, 0.69, 1.13, 1.11, 1.85, 1.13, 0.96, 2.85, 1.39, 1.76, 0.96, 1.39, 1.94, 1.54, 1.54, 1.54, 0.73, 0.73, 1.76, 0.73, 1.39, 0.75, 1.63, 1.76, 0.96]
            work_contacts_sdlog = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.2, 1.2, 1.2, 0.21, 0.21, 1.42, 1.42, 1.42, 0.99, 0.99, 0.99, 0.57, 0.81, 1.29, 1.31, 0.81, 0.99, 1.06, 1.43, 1.42, 0.99, 1.43, 1.56, 1.12, 1.12, 1.12, 0.72, 0.72, 1.42, 0.72, 1.43, 1.4, 1.59, 1.42, 0.99]

            workplace_dynamic_degree_distribution = Array{Distribution,1}(undef, length(work_contacts_meanlog))

            for ii = 1:length(workplace_dynamic_degree_distribution)
                workplace_dynamic_degree_distribution[ii] = Distributions.LogNormal(work_contacts_meanlog[ii],work_contacts_sdlog[ii])
            end
        end

        return workplace_dynamic_degree_distribution
    end
end

"""
Transmission related fns
"""

function transmit_over!(transmission_risk::Float64,
                    timelat::Array{Int64,1},
                    infected_by::Array{Int64,1},
                    output::sim_outputs,
                    states::node_states,
                    probasymp::Float64,
                    rng::MersenneTwister,
                    time::Int64,
                    count::Int64,
                    intervention_set_itr::Int64;
                    infecting_by::Int64=0,
                    contacts_to_check::Array{Int64,1}=Array{Int64,1}(),
                    dynamic_contact::Int64=0,
                    inisol::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    atwork::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    transmission_setting::String = "undefined",
                    social_contact_scaling::Float64 = -1.,
                    random_contact_scaling::Float64 = -1.)
    # Inputs:
    # transmission_risk::Float64 - Probability of transmission given contact made
    # timelat::Array{Int64,1} - Time spent in latent state per node. If zero, then node is susceptible.
    # infected_by::Array{Int64,1} - Record of ID of infectors of each node
    # output::sim_outputs - record of outputs from the simulations
    # states::node_states - status of each node
    # probasymp::Float64 - Asymptomatic probability
    # rng::MersenneTwister - The random number generator
    # time::Int64 - Current timestep
    # count::Int64 - Replicate ID
    # infecting_by::Int64 - Node ID of node transmitting infection
    # contacts_to_check::Int64 - list of nodes IDs to receive infection
    # dynamic_contact::Int64 - Flag for if this infection is over a dynamic contact link
    # inisol::Array{Int64,2} - Flag if nodes are in isolation
    # atwork::Array{Int64,2} - Flag if nodes are at work
    # transmission_setting::String - Specify where potential tranmission event took place. Used to increment output.transmission_setting.

    for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
        infecting_to = contacts_to_check[contact_itr]

        contact_active_check = rand(rng)
        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation or at work if given the arguments inisol / atwork
        if (timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(atwork) == false || atwork[time,infecting_to] == true) &&
            ((social_contact_scaling < 0) || (contact_active_check <= social_contact_scaling)) &&
            ((random_contact_scaling < 0) || (contact_active_check <= random_contact_scaling))

            if rand(rng) < transmission_risk
                timelat[infecting_to] = 1
                infected_by[infecting_to] = infecting_by
                output.num_infected[infecting_by,count,intervention_set_itr] += 1
                states.acquired_infection[infecting_to] = time

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count,intervention_set_itr] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count,intervention_set_itr] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count,intervention_set_itr] += 1

                if transmission_setting == "work"
                    if dynamic_contact==1
                        # Record transmission occurring in work-dynamic setting
                        output.transmission_setting[time+1,count,intervention_set_itr,3] += 1
                        # if this was a dynamic contact, update counter
                        output.dynamic_infection_count[time+1,count,intervention_set_itr] += 1
                    else
                        # Record transmission occurring in work-static setting
                        output.transmission_setting[time+1,count,intervention_set_itr,2] += 1
                    end
                elseif transmission_setting == "social"
                    # Record transmission occurring in social setting
                    output.transmission_setting[time+1,count,intervention_set_itr,1] += 1
                elseif transmission_setting == "household"
                    # Record transmission occurring in household setting
                    output.transmission_setting[time+1,count,intervention_set_itr,4] += 1
                elseif transmission_setting == "other"
                    # Record transmission occurring in other setting
                    output.transmission_setting[time+1,count,intervention_set_itr,5] += 1
                else
                    error("Invalid transmission setting!")
                end
            end
        end
    end
end

function transmit_over_other_workplace!(transmission_risk::Float64,
                                            infected_by::Array{Int64,1},
                                            output::sim_outputs,
                                            states::node_states,
                                            probasymp::Float64,
                                            rng::MersenneTwister,
                                            time::Int64,
                                            count::Int64,
                                            intervention_set_itr::Int64,
                                            infecting_by::Int64,
                                            contacts_to_check::Array{Int64,1},
                                            inisol::Array{Int64,2},
                                            atwork::Array{Int64,2},
                                            network_parameters::network_params)
    # Inputs:
    # transmission_risk::Float64 - Probability of transmission given contact made
    # infected_by::Array{Int64,1} - Record of ID of infectors of each node
    # output::sim_outputs - record of outputs from the simulations
    # states::node_states - status of each node
    # probasymp::Float64 - Asymptomatic probability
    # rng::MersenneTwister - The random number generator
    # time::Int64 - Current timestep
    # count::Int64 - Replicate ID
    # infecting_by::Int64 - Node ID of node transmitting infection
    # contacts_to_check::Int64 - list of nodes IDs to receive infection
    # inisol::Array{Int64,2} - Flag if nodes are in isolation
    # atwork::Array{Int64,2} - Flag if nodes are at work
    # network_parameters::network_params - Quantities to construct the contacts & stores the node properties

    for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
        infecting_to = contacts_to_check[contact_itr]

        # Get information on the contacts workplace
        contact_sector_ID = network_parameters.worker_nodes[infecting_to].sector_ID
        contact_workplace_ID = network_parameters.worker_nodes[infecting_to].workplace_ID
        contact_workplace_info = network_parameters.workplace_info[contact_sector_ID][contact_workplace_ID]

        # Get contacts workplace CS status
        contact_workplace_CS_bool = contact_workplace_info.covid_secure

        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation at at work if given the arguments inisol / atwork
        if (states.timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(atwork) == false || atwork[time,infecting_to] == true) &&
            (contact_workplace_CS_bool==false)

            if rand(rng) < transmission_risk
                states.timelat[infecting_to] = 1
                infected_by[infecting_to] = infecting_by
                output.num_infected[infecting_by,count,intervention_set_itr] += 1
                states.acquired_infection[infecting_to] = time

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count,intervention_set_itr] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count,intervention_set_itr] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count,intervention_set_itr] += 1

                # Record transmission occurring in work-static setting
                output.transmission_setting[time+1,count,intervention_set_itr,2] += 1

            end
        end
    end
end

"""
Functions to set up transmission rates within household, workplace and socially for each individual
"""
function assign_household_transmit_onegroup!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean != 1)
        error("Should only be a single household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 1)
        error("Should only be a single household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    mean_val = transrisk_household_group_mean[1]
    sd_val = transrisk_household_group_sd[1]
    norm_dist = Normal(mean_val,sd_val)

    # Iterate over each individual
    for node_itr = 1:n_nodes
        worker_nodes[node_itr].transrisk_household = rand(rng,norm_dist)
    end

    return nothing
end

# Household risk based on household size
function assign_household_transmit_household_size!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if there is not four groups
    if (n_transrisk_household_group_mean != 4)
        error("Should be four group SAR mean estimates, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 4)
        error("Should be four group SAR standard deviations, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:n_nodes

        # Check number of people in the household
        hh_size = household_contacts_per_node[node_itr] + 1
        if hh_size == 1
            # Sole person in household. No household contacts, set transmission risk to zero.
            worker_nodes[node_itr].transrisk_household = 0
        elseif hh_size == 2
            # Household size two
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[1])
        elseif hh_size == 3
            # Household size three
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[2])
        elseif hh_size == 4
            # Household size four
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[3])
        else
            # Household size five or more
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[4])
        end
    end

    return nothing
end

function assign_household_transmit_multigrouptest!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_household_group = length(transrisk_household_group)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean != 1)
        error("Should only be a single household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 1)
        error("Should only be a single household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end


    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:n_nodes
        if household_contacts_per_node[node_itr] == 0
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[1])
        elseif household_contacts_per_node[node_itr] == 1
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[2])
        else
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[3])
        end
    end

    return nothing
end

function assign_workplace_static_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    workertypes::Int64,
                                    transrisk_static_work_mean::Array{Float64,1},
                                    transrisk_static_work_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_mean = length(transrisk_static_work_mean)
    n_transrisk_sd = length(transrisk_static_work_sd)

    # Throw error if dimensions don't match
    if (n_transrisk_mean != workertypes) ||
        (n_transrisk_sd != workertypes)
        error("Dimension mismatch in transrisk parameter arrays. Please rectify.")
    end

    # Iterate over each individual
    for node_itr = 1:n_nodes
        worker = worker_nodes[node_itr]
        sector = worker.sector_ID
        mean_val = transrisk_static_work_mean[sector]
        sd_val = transrisk_static_work_sd[sector]
        worker_nodes[node_itr].transrisk_static_work = rand(rng,Normal(mean_val,sd_val))
    end

    return nothing
end

function assign_workplace_dynamic_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    workertypes::Int64,
                                    transrisk_dynamic_work_mean::Array{Float64,1},
                                    transrisk_dynamic_work_sd::Array{Float64,1})
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Get number of household groups in use
    n_transrisk_mean = length(transrisk_dynamic_work_mean)
    n_transrisk_sd = length(transrisk_dynamic_work_sd)

    # Throw error if dimensions don't match
    if (n_transrisk_mean != workertypes) ||
        (n_transrisk_sd != workertypes)
        error("Dimension mismatch in transrisk parameter arrays. Please rectify.")
    end

    # Iterate over each individual
    for node_itr = 1:n_nodes
        worker = worker_nodes[node_itr]
        sector = worker.sector_ID
        mean_val = transrisk_dynamic_work_mean[sector]
        sd_val = transrisk_dynamic_work_sd[sector]
        worker_nodes[node_itr].transrisk_dynamic_work = rand(rng,Normal(mean_val,sd_val))
    end

    return nothing
end


function assign_social_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    transrisk_social_mean::Float64,
                                    transrisk_social_sd::Float64)
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Iterate over each individual
    for node_itr = 1:n_nodes
        worker_nodes[node_itr].transrisk_social = rand(rng,Normal(transrisk_social_mean,transrisk_social_sd))
    end

    return nothing
end

function assign_random_transmit!(RNGseed::Int64,
                                    network_parameters::network_params,
                                    transrisk_random_mean::Float64,
                                    transrisk_random_sd::Float64)
# Inputs:
# RNGseed::Int64 - Set the number to seed the random number generator
# network_parameters::network_params - Quantities to construct the contacts & stores the node properties
# household_contacts_per_node::Array{Int64,1} - As described
# transrisk_household_group_mean/sd::Array{Float64,1} -  probability of transmission within a household, can differ per household group

    # Unpack parameters
    @unpack n_nodes, worker_nodes = network_parameters

    # Set the random number generator
    rng = MersenneTwister(RNGseed)

    # Iterate over each individual
    for node_itr = 1:n_nodes
        worker_nodes[node_itr].transrisk_random = rand(rng,Normal(transrisk_random_mean,transrisk_random_sd))
    end

    return nothing
end



"""
Functions to reinitialise states at start of each run
"""
# Node states, household inf delay & CT vars
function reinitialise_node_states!(states::node_states)
    lmul!(0,states.timelat)
    lmul!(0,states.timeinf)
    lmul!(0,states.timesymp)
    lmul!(0,states.asymp)
    lmul!(0,states.lattime)
    lmul!(0,states.timeisol)
    lmul!(0,states.symp_timeisol)
    lmul!(0,states.timeisol_CTcause)
    lmul!(0,states.hh_isolation)
    lmul!(0,states.delay_adherence)
    lmul!(0,states.acquired_infection)
end

function reinitialise_daily_record_arrays!(contacts::contacts_struct)
    lmul!(0,contacts.daily_record_atworkplace)
    lmul!(0,contacts.daily_record_inisol)
end

function reinitialise_workplace_params!(workplace_info::Array{Array{workplace_params,1},1})

    # Get number of sectors in use
    n_sectors = length(workplace_info)

    # Iterate over each sector.
    # Within each sector, iterate over each workplace
    # Update field values
    for sector_itr = 1:n_sectors
        n_workplaces = length(workplace_info[sector_itr])
        for workplace_itr = 1:n_workplaces
            #workplace_info[sector_itr][workplace_itr].covid_secure = false
            workplace_info[sector_itr][workplace_itr].workplace_open = true
        end
    end

    return nothing
end

function reinitialise_CT_vars!(CT_vars::contact_tracing_vars,cmax::Int64, rng::MersenneTwister,
    CT_parameters::CT_params, delay_adherence::Array{Int64,1},
    csum_test_result_delay::Array{Float64,1},max_test_result_delay::Int64)

@unpack CT_days_before_symptom_included, CT_engagement = CT_parameters

    # Reset vector tracking symptomatic cases (positive confirmed or untested)
    lmul!(0,CT_vars.Symp_cases_per_household_pos_or_unknown)

    # Reset vector tracking the latest isolation release time due to household member cases
    lmul!(0,CT_vars.hh_isolation_release_time)

    # Variables for waiting for test results
    lmul!(0,CT_vars.Time_to_test_result)
    CT_vars.Time_to_test_result .-= 1 # Reset so all values are -1

    # Repopulate Boolean vector stating whether a false negative test result would be returned
    # and the number of days relevant for contact tracing
    lmul!(0,CT_vars.relevant_prev_days_for_CT)
    for ii = 1:cmax

      # For each worker, initialise CT_vars.Test_result_false_negative as false
      CT_vars.Test_result_false_negative[ii] = false

      # For each worker, check if they engage with contact tracing
      engage_with_CT_rand = rand(rng)
      if engage_with_CT_rand < CT_engagement # engage with contact tracing
          CT_vars.Engage_with_CT[ii] = true
      else # do not engage with contact tracing
          CT_vars.Engage_with_CT[ii] = false
      end

      # Get amount of days to be looked back over
      # Have upper bound of 7 days post symp
      # if we put in reporting delay, needs to be above this
      CT_vars.relevant_prev_days_for_CT[ii] = min(CT_days_before_symptom_included + delay_adherence[ii],
                                          CT_days_before_symptom_included + 7)
    end

    # Repopulate time until test result received for each individual
    lmul!(0,CT_vars.CT_delay_until_test_result)
    for node_itr = 1:cmax
                      CT_vars.CT_delay_until_test_result[node_itr] = draw_sample_from_pmf(csum_test_result_delay,
                                                                                rng;
                                                                                idx_offset = 1)
    end

    # Set up vector of vectors for storing IDs of those to be contacted in CT
    CT_vars.Inds_to_be_contacted = Array{Array{Int64,1},1}(undef,cmax)

    # Initialise array to keep track of whether an infected recalls their infector
    lmul!(0,CT_vars.Recall_infector)
end

"""
Misc. fns
"""

function draw_sample_from_pmf(csum_pmf::Array{Float64,1},
                                rng::MersenneTwister;
                                idx_offset::Int64 = 0)
# Inputs:
# val_to_update::Int64 - Entry sampled value will be assigned to
# csum_pmf::Array{Float64,1} - Cumulative summed probability mass function. Used to draw value from.
# rng::MersenneTwister - The random number generator
# idx_offset::Int64 = 0 - Links bin index to the quantity value

    # Get number of elements in the pmf
    n_bins = length(csum_pmf)

    # Initialise output value
    val_to_update = 0

    # Draw random number
    # Set delay in adherence/symptoms becoming known to household
    # Find interval random number resides in
    r = rand(rng)
    allocated_flag = false # Intialise allocation flag. Switch to true when allocation done.
    bin_idx = 1   # Current interval being checked
    while (allocated_flag == false)
        if r <= csum_pmf[bin_idx]
            # Assign selected value
            # Subtract idx_offset
            val_to_update = bin_idx - idx_offset

            # Update allocation flag
            allocated_flag = true
        else
            # r does not reside in this interval. Update bin index.
            bin_idx += 1

            # Error check, if not assigned value after checked final bin value
            if bin_idx > n_bins
                error("bin_idx is now $bin_idx. The pmf only has $n_bins bins. Terminating programme.")
            end
        end
    end

    return val_to_update::Int64
end

function set_infection_related_times!(time_to_symps::Array{Int64,1},states::node_states,
    isolation::Int64,adherence::Float64,csum_delay_adherence::Array{Float64,1},
    d_incub::Distribution,cmax::Int64,rng::MersenneTwister)

    time_to_symps .= ceil.(rand(rng,d_incub,cmax)) # time to symptoms
    # (for asymptomatics, the same from a silent start of "symptoms")

    # iterate over nodes to set lattime and hh_isolation
    for node_itr = 1:cmax
        # lattime is the time from infection to infectiousness
        if time_to_symps[node_itr]-states.inftime<1
            states.lattime[node_itr] = 1  # Infectiousness can begin the day after becoming infected
        else
            states.lattime[node_itr] = time_to_symps[node_itr]-states.inftime
        end


        if isolation==1
            p1 = rand(rng)
            if p1 < adherence # those who adhere will isolate when they get symptoms
                states.hh_isolation[node_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
            end

            # Draw random number
            # Set delay in adherence/symptoms becoming known to household
            # Find interval random number resides in
            states.delay_adherence[node_itr] = draw_sample_from_pmf(csum_delay_adherence,
                                                                    rng;
                                                                    idx_offset = 1)
        end

    end
end

# function increment_infection_process!(states::node_states,
#     output::sim_outputs, worker_nodes::Array{worker_params,1},
#     workplace_closure_active,
#     infected_by,
#     probasymp::Float64,
#     rng::MersenneTwister,
#     time::Int64,
#     count::Int64,
#     household_contacts::Array{Array{Int64,1},1},
#     household_contacts_per_node::Array{Int64,1},
#     contact_tracing_active,
#     CT_vars::contact_tracing_vars,
#     workplace_memory,
#     temp::temp_struct,
#     undefined_array::Array{Int64,2})
#
#     # If come to the end of latent time, move to infectious time etc
#     for node_itr = 1:cmax
#         # if the node has reached the end of latent infection
#         if states.timelat[node_itr]>states.lattime[node_itr]
#             # move to being infectious
#             states.timelat[node_itr] = -1
#             states.timeinf[node_itr] = 1
#
#             # Increment time series counts
#             output.numinf[time,count] += 1
#             output.newinf[time,count] += 1
#
#             # check if new infected will be asymptomatic
#             if states.asymp[node_itr] > 0
#                 output.newasymp[time,count] += 1
#             end
#
#             # check if it is worker that is newly infected
#             if worker_nodes[node_itr].returned_to_work==1
#                 output.workersinf[time,count] += 1
#
#                 # Count workers that are asymptomatically infected
#                 if states.asymp[node_itr] > 0
#                     output.workersasymp[time,count] += 1
#                 elseif workplace_closure_active==true
#                     # if the newly infected worker is symptomatic, add to
#                     # the workplace memory
#                     temp.workertype_ID = worker_nodes[node_itr].sector_ID
#                     temp.workplace_ID = worker_nodes[node_itr].workplace_ID
#                     workplace_memory[temp.workertype_ID][temp.workplace_ID,WP_memory_slot] += 1
#                 end
#             end
#
#             # infect rest of household
#             # transmit over household_contacts[node_itr]
#             # checking that contacts are susceptible
#             transmit_over!(1.,states.timelat,infected_by,output,probasymp,
#                     states.asymp,rng,time,count,
#                     infecting_by=node_itr,
#                     contacts_to_check=household_contacts[node_itr],
#                     inisol=undefined_array,
#                     atwork=undefined_array)
#         end
    #     # Update node disease state time vectors
    #     if states.timeinf[node_itr]>states.inftime
    #         # the node becomes symptomatic (if they develop symptoms)
    #         states.timeinf[node_itr] = -1
    #         states.timesymp[node_itr] = 1
    #
    #         # Increment time series counts
    #         output.numrep[time,count] += 1
    #
    #         # Check if index case are symptomatic & would have zero adherence delay
    #         if (states.asymp[node_itr] == 0) && (states.delay_adherence[node_itr]==0)
    #             # Check if infected will isolate
    #             if (states.hh_isolation[node_itr]==1)
    #                 states.symp_timeisol[node_itr] = 1
    #
    #                 # Set that the unit has reported infection this timestep
    #                 states.rep_inf_this_timestep[node_itr] = 1
    #             end
    #
    #             # Irrespective of whether index case self-isolates,
    #             # adherent members of their household may also isolate.
    #             for hh = 1:household_contacts_per_node[node_itr]
    #                 temp.contact_ID = household_contacts[node_itr][hh]
    #                 if (states.hh_isolation[temp.contact_ID]==1) && (states.symp_timeisol[temp.contact_ID]==0) # Individual not already symptomatic themselves
    #                     states.timeisol[temp.contact_ID] = 1
    #                 end
    #             end
    #         end
    #
    #         # If contact tracing active, increase number of symptomatic infections
    #         # in household by one
    #         if contact_tracing_active == true
    #             if (states.asymp[node_itr] == 0) # Check case is symptomatic
    #                 current_node_household_ID = worker_nodes[node_itr].household_ID
    #                 CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] += 1
    #             end
    #         end
    #     end
    #
    #     # Check if node, if having a delayed adherence, begins adherence on current day
    #     if (states.timesymp[node_itr] > 1)&&((states.timesymp[node_itr]-1)==states.delay_adherence[node_itr]) # Condition for node beginning adherence on current day & has been symptomatic for at least one day
    #         if states.asymp[node_itr] == 0 # Check node is symptomatic and will adhere
    #             if states.hh_isolation[node_itr]==1 # Check node will adhere
    #                 states.symp_timeisol[node_itr] = 1 + states.delay_adherence[node_itr]
    #
    #                 # Set that the unit has reported infection this timestep
    #                 states.rep_inf_this_timestep[node_itr] = 1
    #             end
    #
    #             # Household members now aware of index case having symptoms.
    #             # Adherent members of their household may also now isolate,
    #             # assuming infected displays symptoms
    #             # Note they are delayed in isolating, in line with delay
    #             # of the index case
    #             for hh = 1:household_contacts_per_node[node_itr]
    #                 temp.contact_ID = household_contacts[node_itr][hh]
    #                 if (states.hh_isolation[temp.contact_ID]==1) && (states.symp_timeisol[temp.contact_ID]==0) # Individual not already symptomatic themselves
    #                         # Individual shortens isolation by length of time since
    #                         # unwell individual began displaying symptoms
    #                     states.timeisol[temp.contact_ID] = 1 + states.delay_adherence[node_itr]
    #                 end
    #             end
    #         end
    #     end
    #
    #     # Check if node has reached end of symptom period
    #     if states.timesymp[node_itr]>states.symptime
    #         states.timesymp[node_itr] = -1
    #
    #         # If contact tracing active and case was symptomatic,
    #         # decrease number of symptomatic infections in household by one
    #         if contact_tracing_active == true
    #             # Check case is symptomatic & not returned a false negative (if false negative, has already been subtracted)
    #             if (states.asymp[node_itr] == 0) && (CT_vars.Test_result_false_negative[node_itr] == false)
    #                 current_node_household_ID = worker_nodes[node_itr].household_ID
    #                 CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1
    #
    #                 # Error check
    #                 if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
    #                     error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
    #                 end
    #
    #             end
    #         end
    #     end
    # end
# end
