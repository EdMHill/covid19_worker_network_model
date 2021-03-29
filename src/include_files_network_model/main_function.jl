"""
Main function
"""

# Run outbreak on network
function worker_pattern_network_run(RNGseed::Int64,
                                        cmax::Int64,
                                        ton::Int64,
                                        toff::Int64,
                                        infection_parameters::infection_params,
                                        sameday::Int64,
                                        seed_initial_states_fn::Function,
                                        countfinal::Int64,
                                        endtime::Int64,
                                        contact_tracing_active::Bool,
                                        CT_parameters::CT_params,
                                        network_parameters::network_params,
                                        workplace_generation_parameters::workplace_generation_params,
                                        workplace_closure_active::Bool,
                                        intervention_list_configs::Array{Array{intervention_struct,1},1},
                                        runset::String;
                                        intervention_fns::Array{Function} = Array{Function}(undef),
                                        assign_household_transrisk_fn::Function = assign_household_transmit_onegroup!,
                                        assign_workplace_static_transrisk_fn::Function = assign_workplace_static_transmit!,
                                        assign_workplace_dynamic_transrisk_fn::Function = assign_workplace_dynamic_transmit!,
                                        assign_social_transrisk_fn::Function = assign_social_transmit!,
                                        assign_random_transrisk_fn::Function = assign_random_transmit!,
                                        include_plot::Bool = false,
                                        include_statistics::Bool = false,
                                        include_dd::Bool = false,
                                        include_dd_dynamic_work::Bool = false,
                                        include_dd_dynamic_social::Bool = false)

# Inputs:
# RNGseed - Sets the random number generator
# cmax::Int64 - Number of nodes
# ton::Int64, toff::Int64 - Work pattern vars. Days on and days off, for sameday==3 ton weeks on followed by toff weeks off
# infection_parameters::infection_params - Parameter structure for infection params
# sameday::Int64 - Flag.  If sameday=0, all workers are at work on the same set of consecutive days. If sameday=1, workers go to work on a random set of consecutive days. If sameday=3, workers go to work on the same number of days, but scattered randomly throughout the week.
# seed_initial_states_fn::Function - Sets amount of nodes to be seeded as latent/asymp/symp/recovered.
# countfinal::Int64, endtime::Int64 - % Simulation replicates & timesteps for each replicate
# endtime::Int64 - Number of timesteps each simulation replicate to be performed for
# contact_tracing_active::Bool - Set if contact tracing is active or not (Bool type variable)
# CT_parameters::CT_params - Parameter structure for contact tracing params
# network_parameters::network_params - Parameter structure for network params
# workplace_generation_parameters::workplace_generation_params - Parameter structure for workplace generation params
# workplace_closure_active::Bool - Whether workplace closures are in operation or not
# intervention_list_configs - Specify use of any additional, trigger interventions
# runset::String - Name of the scenario being run. Used to bring interventions into effect for specified scenarios.
# assign_household_transrisk_fn - Specify assignment of individuals to household groups with differing household transmission risk

# For parameters within the parameter structures, see "include_files_network_model/parametertypes.jl"

##  OUTLINE OF THE CODE STRUCTURE
# - Unpack required variables
# - Set the random number generator
# - Initialise workplaces and workers
# - Generate contacts within workplaces
# - Generate social contacts (workdays and non-workdays)
# - Generate dynamic contacts
# - Generate household specific transmission
# - Initialise variables
# - Iterate over replicates
#        -- Set course of infection times
#        -- Reset contact tracing variables
#        -- Iterate over time
#               --- Clear workplace memory for this timestep
#               --- Assign outputs
#               --- Increment counters
#               --- Increment infection process (if come to the end of latent time, move to infectious time etc.)
#                       -> Includes household infection loop
#               --- Transmit infections
#               --- Perform contact tracing
#               --- Assign isolation outputs
#               --- Close workplaces
#               --- Run interventions

##

     """
     Unpack required variables
     """
    # if contact_tracing_active==true
        @unpack CT_engagement, CT_delay_until_test_result_pmf, CT_days_before_symptom_included, test_detection_prob_vec,
            CT_caused_isol_limit, dynamic_contacts_recalled_propn, social_contacts_recalled_propn, prob_backwards_CT,
            perform_CT_from_infector, infector_engage_with_CT_prob, contact_tracing_active, workplace_closure_active = CT_parameters
    # end
    # if workplace_closure_active==true
        @unpack workplace_CT_memory, workplace_CT_threshold, time_WC = CT_parameters
    # end
    @unpack asymp_trans_scaling_dist, iso_trans_scaling, CS_scale_transrisk,
            transrisk_household_group_mean, transrisk_household_group_sd,
            transrisk_static_work_mean, transrisk_static_work_sd,
            transrisk_dynamic_work_mean, transrisk_dynamic_work_sd,
            transrisk_social_mean, transrisk_social_sd,
            transrisk_random_mean, transrisk_random_sd,
            probasymp_dist, isolation, symp_isoltime, household_isoltime, adherence,
            n_social_mean_workday,n_social_mean_nonworkday,
            d_incub, dist_infectivity, delay_adherence_pmf,
            delay_household_infection_pmf, recov_propn = infection_parameters
    @unpack dynamic_conts_mean, dynamic_conts_sd, CS_active_flag = network_parameters
    @unpack workertypes = workplace_generation_parameters

    """
    Set the random number generator
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise workplaces and workers
    """

    # generate_workplaces_and_allocate_workers in network_generation_fns.jl
    # households [nodes x workers x properties]: dim1--nodes, dim2--array containing workers, dim3--tuple of worker properties
    # worker properties: (0/1: returned to work (1) or not (0), worker group id, id of workplace within worker group)
    # workplace_sizes [workertype x workplace]: sizes of workplaces within each worker group
    # workplace_info: Parameter structure to store details on each workplace
    @time worker_nodes::Array{worker_params,1},
    workplace_sizes::Array{Array{Int64,1},1},
    workplace_info::Array{Array{workplace_params,1},1},
    nodes_by_workplace::Array{Array{Array{Int64,1},1},1} = generate_workplaces_and_allocate_workers(cmax, workplace_generation_parameters,RNGseed,rng,CS_active_flag)
    network_parameters.worker_nodes = worker_nodes
    network_parameters.workplace_sizes = workplace_sizes
    network_parameters.workplace_info = workplace_info

    workplace_size_1D = collect(Iterators.flatten(workplace_sizes))
    max_internal_workplace_contacts = sum(workplace_size_1D.*(workplace_size_1D.-1))/2
    println("Max internal workplace contacts: $max_internal_workplace_contacts")

    # Initialise sector open/closure status
    sector_open = Array{Bool,1}(undef,workertypes)
    for sector_itr = 1:workertypes
        sector_open[sector_itr] = true
    end
    network_parameters.sector_open = sector_open

    println("cmax inside: $cmax")

    """
    Generate contacts within workplaces & households
    """
    # Generate contact network within workplaces, in network_generation_fns.jl
    @time contacts::contacts_struct,
        social_contacts::Array{Array{Int64,1},1},
        household_contacts::Array{Array{Int64,1},1},
        work_contacts_same_workplace_per_node::Array{Int64,1},
        work_contacts_other_workplace_per_node::Array{Int64,1},
        household_contacts_per_node::Array{Int64,1},
        n_households::Int64 = generate_contacts(cmax,
                                                endtime,
                                                network_parameters,
                                                workplace_generation_parameters,
                                                nodes_by_workplace,
                                                RNGseed,rng)

    println("worker_nodes: $(worker_nodes[1:10])")
    println("workplace_sizes: $(workplace_sizes)")
    println("worker_nodes.household_ID: $(worker_nodes[1].household_ID)")
    println("n_households: $(n_households)")


    if include_plot
        fig = plot_network(worker_nodes,
                        contacts,
                        social_contacts,
                        household_contacts)
        return fig
    end

    """
    Generate social contacts (workdays and non-workdays)
    """
    # generate_social_contacts_each_day in network_generation_fns.jl
    # Set up social contacts made on work days for each given day
    # & set up social contacts made on non-work days for each given day
    # workday_social_contacts_by_day,

    if network_parameters.network_generation_method_dynamic_social == "ER"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                    RNGseed,
                                                                                    cmax,
                                                                                    endtime,
                                                                                    social_contacts,
                                                                                    contacts.social_contacts_per_node,
                                                                                    n_social_mean_workday,
                                                                                    n_social_mean_nonworkday)
    elseif network_parameters.network_generation_method_dynamic_social == "configuration"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                RNGseed,
                                                                                cmax,
                                                                                endtime,
                                                                                social_contacts,
                                                                                contacts.social_contacts_per_node,
                                                                                network_parameters.social_workday_dd,
                                                                                network_parameters.social_nonworkday_dd,
                                                                                network_parameters.cluster_social_contacts)

    elseif network_parameters.network_generation_method_dynamic_social == "cluster"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                RNGseed,
                                                                                cmax,
                                                                                endtime,
                                                                                social_contacts,
                                                                                network_parameters.social_workday_dd,
                                                                                network_parameters.social_nonworkday_dd)

    elseif network_parameters.network_generation_method_dynamic_social == "repeated"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                RNGseed,
                                                                                cmax,
                                                                                endtime,
                                                                                social_contacts,
                                                                                network_parameters.social_workday_dd,
                                                                                network_parameters.group_limit,
                                                                                network_parameters.dynamic_time_frame)

    elseif network_parameters.network_generation_method_dynamic_social == "groups"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                RNGseed,
                                                                                cmax,
                                                                                endtime,
                                                                                network_parameters.n_groups_per_day_distribution,
                                                                                network_parameters.group_limit,
                                                                                network_parameters.dynamic_time_frame)
        contacts.social_contacts_per_node .+= 1

    elseif network_parameters.network_generation_method_dynamic_social == "fixed_daily_degree"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day(rng,
                                                                                RNGseed,
                                                                                cmax,
                                                                                endtime,
                                                                                12,
                                                                                network_parameters.group_limit,
                                                                                network_parameters.dynamic_time_frame)
        contacts.social_contacts_per_node .+= 1
    else
        println("Invalid network generation method!")
    end

    contacts.workday_social_contacts_by_day = workday_social_contacts_by_day
    contacts.nonworkday_social_contacts_by_day = nonworkday_social_contacts_by_day

    """
    Generate dynamic contacts
    """
    # Set up dynamic worker contacts, in network_generation_fns.jl

    if network_parameters.network_generation_method == "ER"
        dynamic_worker_contacts = generate_dynamic_worker_contacts(RNGseed,
                                                                    cmax,
                                                                    endtime,
                                                                    worker_nodes,
                                                                    dynamic_conts_mean,
                                                                    dynamic_conts_sd)

    elseif network_parameters.network_generation_method == "configuration"
        dynamic_worker_contacts = generate_dynamic_worker_contacts(RNGseed,
                                                                cmax,
                                                                endtime,
                                                                worker_nodes,
                                                                network_parameters.workplace_dynamic_degree_distribution,
                                                                network_parameters.max_contacts_work_dynamic)
    else
        println("Invalid network generation method!")
    end

    contacts.dynamic_worker_contacts = dynamic_worker_contacts

    """
    Generate random dynamic contacts
    """
    dynamic_random_contacts = generate_random_contacts(RNGseed,
                                                        cmax,
                                                        endtime,
                                                        network_parameters.prob_random_contact)

    contacts.dynamic_random_contacts = dynamic_random_contacts


    if include_statistics
        @time gcc_static = calc_static_gcc(worker_nodes,
                        contacts,
                        social_contacts,
                        household_contacts)

        states = node_states(cmax=cmax)
        @time populate_atwork!(states.atwork,sameday,ton,toff,rng)

        @time gcc_dynamic = calc_dynamic_gcc(worker_nodes,
                        contacts,
                        social_contacts,
                        household_contacts,
                        states.atwork,
                        14)

        return gcc_static, gcc_dynamic
    end

    if include_dd
        @time dd_work, dd_social = calc_dd_static(worker_nodes,
                                            contacts,
                                            social_contacts,
                                            household_contacts,
                                            workertypes)

        return dd_work, dd_social
    end

    if include_dd_dynamic_work
        @time dd_dynamic_work = calc_dd_dynamic_work(worker_nodes,
                                                    contacts,
                                                    workertypes,
                                                    endtime)

        return dd_dynamic_work
    end

    if include_dd_dynamic_social
        @time dd_dynamic_social, proportion_friends_met = calc_dd_dynamic_social(worker_nodes,
                                                    contacts,
                                                    social_contacts,
                                                    endtime)

        return dd_dynamic_social, proportion_friends_met
    end


    """
    Generate transmission risks
    """
    # Relevant functions are listed in "include_files_network_model/additional_fns.jl"
    assign_household_transrisk_fn(RNGseed,
                                        network_parameters,
                                        household_contacts_per_node,
                                        transrisk_household_group_mean,
                                        transrisk_household_group_sd)

    assign_workplace_static_transrisk_fn(RNGseed,
                                        network_parameters,
                                        workertypes,
                                        transrisk_static_work_mean,
                                        transrisk_static_work_sd)

    assign_workplace_dynamic_transrisk_fn(RNGseed,
                                        network_parameters,
                                        workertypes,
                                        transrisk_dynamic_work_mean,
                                        transrisk_dynamic_work_sd)

    assign_social_transrisk_fn(RNGseed,
                                network_parameters,
                                transrisk_social_mean,
                                transrisk_social_sd)

    assign_random_transrisk_fn(RNGseed,
                                network_parameters,
                                transrisk_random_mean,
                                transrisk_random_sd)

    """
    Initialise storage arrays
    """

    # Initialise output structure
    output = sim_outputs(endtime=endtime,countfinal=countfinal,cmax=cmax,
                            n_intervention_sets = length(intervention_list_configs))

    # Initialise vectors used in each replicate
    dayon = zeros(Int64,cmax)
    ccount = zeros(Int64,cmax)
    pflag = zeros(Int64,cmax)

    # Initialise node states
    states = node_states(cmax=cmax)
    time_to_symps = zeros(Int64,cmax)
    infected_by = zeros(Int64,cmax)

    # delay_adherence = zeros(Int64,cmax)     # Individial may not report symptoms immediately.
    max_delay_adherence = length(delay_adherence_pmf) # Vars to be used when allocating reporting delay
    csum_delay_adherence = cumsum(delay_adherence_pmf)

    # Initialise undefined array
    # Used transmit_over! function
    undefined_array = Array{Int64,2}(undef,0,0)

    # Initialise variables used to allocate household infection delay
    # & perform error check ensuring the pmf sums to 1
    if (sum(delay_household_infection_pmf) > 1+eps()) || (sum(delay_household_infection_pmf) < 1-eps())
        error("delay_household_infection_pmf does not sum to 1. Check values.")
    end
    csum_delay_household_inf = cumsum(delay_household_infection_pmf)

    """
    If required, initialise contract tracing related variables
    """
    # If contact tracing in use, create variables
    # if contact_tracing_active == true

        CT_vars = contact_tracing_vars(cmax=cmax,endtime=endtime,n_households=n_households)

        max_test_result_delay = length(CT_delay_until_test_result_pmf) # Vars to be used when allocating test delay times
        csum_test_result_delay = cumsum(CT_delay_until_test_result_pmf)

        # Array to track amount of ppl isolating as a result of contact tracing.
        # Row per timestep, column per replicate.
        num_isolating_CTcause = zeros(endtime,countfinal)
    # else
    #     # otherwise we'll make a zero version of the variable
    #     CT_vars = contact_tracing_vars()
    # end

    """
    If required, initialise workplace closure variables
    """
    workplace_memory = Array{Array{Int64,2},1}(undef, workertypes) # initialise the memory for each workplace

    # if workplace_closure_active
        work_closed_time = Array{Array{Int64,1},1}(undef, workertypes) # initialise the timer for work closure
        workplace_thresholds = Array{Array{Int64,1},1}(undef, workertypes) # initialise the threshold for work closure
        # workplace_memory = Array{Array{Int64,2},1}(undef, workertypes) # initialise the memory for each workplace
        num_workplaces = zeros(Int64,workertypes)
        for worktypeID = 1:workertypes
            num_workplaces[worktypeID] = length(workplace_sizes[worktypeID])
            work_closed_time[worktypeID] = zeros(num_workplaces[worktypeID])
            workplace_memory[worktypeID] = zeros(num_workplaces[worktypeID],workplace_CT_memory)
            workplace_thresholds[worktypeID] = zeros(num_workplaces[worktypeID])
            for work_ID = 1:num_workplaces[worktypeID]
                workplace_thresholds[worktypeID][work_ID] = ceil(Int64,workplace_sizes[worktypeID][work_ID]*workplace_CT_threshold)
            end
        end
    # end

    """
    If required, initialise triggered intervention variables
    """
    # Check if any interventions were specified
    if isassigned(intervention_fns)
        # Number of intervention sets provided is number of rows of intervention_fns
        n_intervention_fns = size(intervention_fns,1)
    end

    # Store preintervention variables
    if length(intervention_list_configs[1]) > 0
        workplace_generation_parameters_preintervention = deepcopy(workplace_generation_parameters)
        states_preintervention = deepcopy(states)
        CT_parameters_preintervention = deepcopy(CT_parameters)
        infection_parameters_preintervention = deepcopy(infection_parameters)
        network_parameters_preintervention = deepcopy(network_parameters)
    end
    """
    Run different intervention sets
    """
    # Find number of intervention sets
    n_intervention_sets = length(intervention_list_configs)

    for intervention_set_itr = 1:n_intervention_sets

        println("Intervention set: $intervention_set_itr")

        intervention_list = intervention_list_configs[intervention_set_itr]

        # Find number and times of specified interventions
        n_interventions = length(intervention_list)
        if n_interventions > 0
            intervention_times = [intervention_list[i].start_time for i=1:n_interventions]
        end

        """
        Run replicates
        """
        # Perform countfinal number of replicates
        for count=1:countfinal

            """
            Set the RNG
            """
            rng = MersenneTwister(RNGseed+count)

            """
            Initialisation phase
            """

            # Re-initiailise node related arrays
            lmul!(0,dayon)
            lmul!(0,ccount)
            lmul!(0,pflag)

            # Reset to pre-intervention conditions
            if n_interventions > 0
                workplace_generation_parameters = deepcopy(workplace_generation_parameters_preintervention)
                states = deepcopy(states_preintervention)
                CT_parameters = deepcopy(CT_parameters_preintervention)
                infection_parameters = deepcopy(infection_parameters_preintervention)
                network_parameters = deepcopy(network_parameters_preintervention)

                @unpack CT_engagement, CT_delay_until_test_result_pmf, CT_days_before_symptom_included, test_detection_prob_vec,
                    CT_caused_isol_limit, dynamic_contacts_recalled_propn, social_contacts_recalled_propn, prob_backwards_CT,
                    perform_CT_from_infector, infector_engage_with_CT_prob, contact_tracing_active, workplace_closure_active = CT_parameters

                @unpack worker_nodes = network_parameters
            end

            # Construct & populate arrays signifiying when at workplace
            populate_atwork!(states.atwork,sameday,ton,toff,rng)

            # Reinitialise time series vectors
            reinitialise_node_states!(states)

            # Reinitialise daily record arrays
            reinitialise_daily_record_arrays!(contacts)

            # Reinitialise workplace_params
            reinitialise_workplace_params!(workplace_info)

            # Reinitialise infected_by & time_to_symps
            lmul!(0,infected_by)
            lmul!(0,time_to_symps)


            """
            Set course of infection times
            """
            # set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence
            set_infection_related_times!(time_to_symps,states,isolation,adherence,csum_delay_adherence,d_incub,cmax,rng)

            """
            Seed non-susceptible disease states
            """
            # Draw asymptomatic probability for current replicate
            probasymp = rand(rng,probasymp_dist)

            # Draw relative infectiousness of an asymptomatic for current replicate
            asymp_trans_scaling = rand(rng,asymp_trans_scaling_dist)

            println("probasymp: $probasymp. asymp_trans_scaling: $asymp_trans_scaling")

            # Sets latent, asymptomatic, symptomatic, recovered nodes
            n_initial_latent::Int64,
            n_initial_asymp::Int64,
            n_initial_symp::Int64,
            n_initial_recovereds::Int64 = seed_initial_states_fn(rng,
                                                            cmax,
                                                            states,
                                                            probasymp,
                                                            infected_by,
                                                            recov_propn)

           """
           Update output time series with initial conditions
           """
           # Update time series for latent & infecteds after assigning initial
           # infecteds
           output.numlat[1,count,intervention_set_itr] = n_initial_latent
           output.numinf[1,count,intervention_set_itr] = n_initial_asymp + n_initial_symp

           # Update prevalences
           output.prevlat[1,count,intervention_set_itr] = n_initial_latent
           output.prevasymp[1,count,intervention_set_itr] = n_initial_asymp
           output.prevsymp[1,count,intervention_set_itr] = n_initial_symp
           output.prevrec[1,count,intervention_set_itr] = n_initial_recovereds

           """
           Reset contact tracing variables
           """
            # If required, set up and/or reinitialise contact tracing related variables
            # if contact_tracing_active == true
                reinitialise_CT_vars!(CT_vars, cmax, rng, CT_parameters, states.delay_adherence,csum_test_result_delay,max_test_result_delay)
            # end

           # Initialise counter for being able to identify the infector of an infectee
           recall_infector_count = 0

           # Counter for identified infector engaging in contact tracing
           # (having not participated in CT before)
           infector_trace_count = [0]

           """
           Generate random number arrays
           """
           r_symp_test = rand(rng,cmax,endtime)
           r_backwards_CT = rand(rng,cmax,endtime)

            """
            Run single replicate
            """
            for time=1:endtime

                # Initial timepoint is for initial conditions
                # Set row to accessed in output arrays for this timestep
                output_time_idx = time + 1

                 """
                 Reinitialise variables at start of timestep
                 """
                # Reinitialise timestep specific values
                lmul!(0,states.rep_inf_this_timestep)

                # reinitalise the current workplace_memory slot
                # if workplace_closure_active == true
                    WP_memory_slot = mod1(time,workplace_CT_memory)
                    for worktypeID = 1:workertypes
                        # Iterate over each workplace for current sector type
                        n_workplaces = num_workplaces[worktypeID]
                        for workplace_itr = 1:n_workplaces
                            workplace_memory[worktypeID][workplace_itr,WP_memory_slot] = 0
                        end
                    end
                # end

                """
                Implement any interventions
                """
                if n_interventions > 0
                    if time ∈ intervention_times
                        intervention_id = findfirst(intervention_times.==time)
                        affect_intervention!(intervention_list[intervention_id],
                                             rng,
                                             RNGseed,
                                            cmax,
                                            ton,
                                            toff,
                                            infection_parameters,
                                            infection_parameters_preintervention,
                                            sameday,
                                            countfinal,
                                            endtime,
                                            CT_parameters,
                                            CT_vars,
                                            contacts,
                                            network_parameters,
                                            nodes_by_workplace,
                                            workplace_generation_parameters,
                                            workplace_closure_active,
                                            states,
                                            workplace_thresholds,
                                            workertypes,
                                            workplace_sizes,
                                            assign_workplace_static_transrisk_fn,
                                           assign_workplace_dynamic_transrisk_fn,
                                           assign_social_transrisk_fn,
                                           assign_random_transrisk_fn)

                        if "contact_tracing" ∈ intervention_list[intervention_id].effects
                            @unpack CT_engagement, CT_delay_until_test_result_pmf, CT_days_before_symptom_included, test_detection_prob_vec,
                                CT_caused_isol_limit, dynamic_contacts_recalled_propn, social_contacts_recalled_propn, prob_backwards_CT,
                                perform_CT_from_infector, infector_engage_with_CT_prob, contact_tracing_active, workplace_closure_active = CT_parameters
                        end
                        if "transrisk" ∈ intervention_list[intervention_id].effects
                            @unpack worker_nodes = network_parameters
                        end
                        if "COVID_secure" ∈ intervention_list[intervention_id].effects
                            @unpack CS_active_flag, workplace_info, worker_nodes = network_parameters
                            @unpack CS_scale_transrisk = infection_parameters
                        end

                        # For specified runsets, with intervention including introduction of isolation guidance
                        # Check symptomatic status of all individuals
                        # Iterate over each symptomatic.
                        # - If they adhere, they isolate
                        # - Check household members. If they adhere (and not symptomatic themselves),
                        #       they enter household isolation tracker.
                        if runset ∈ ["synchronised_changedays_intervention","variable_changedays_intervention",
                                    "synchronised_changedays_intervention_full_adherence","variable_changedays_intervention_full_adherence",
                                    "synchronised_changedays_intervention_low_adherence","variable_changedays_intervention_low_adherence",
                                    "workpercent_intervention","workpercent_intervention_full_adherence","workpercent_intervention_low_adherence",
                                    "CS_intervention","CS_intervention_no_isol","CS_intervention_full_isol",
                                    "adherence_intervention"]
                            for node_itr = 1:cmax
                                # Check symptomatic status of all individuals
                                # Two conditions: Not asymptomatic AND is in symptomatic phase of infection (if symptoms displayed)
                                if (states.asymp[node_itr] == 0) && (states.timesymp[node_itr] > 0)
                                    # For each symptomatic.
                                    # - If the person adheres to new guidance, they now enter symptomatic isolation
                                    # - Check household members. If they adhere (and not symptomatic themselves),
                                    #       they enter household isolation tracker.

                                    # Index case isolation check
                                    if states.hh_isolation[node_itr] == 1

                                        # Isolation due to presence of symptoms
                                        # Set tracker to match time elapsed since symptom onset
                                        time_left_in_symp_isol = (symp_isoltime - states.timesymp[node_itr])
                                            # No +1, as increment on states.timesymp[node_itr] occurs after this loop
                                        for time_itr = 1:time_left_in_symp_isol
                                            array_time_idx = time + time_itr
                                            if array_time_idx <= (endtime + 1)
                                                    states.symp_isolation_array[node_itr,array_time_idx] = 1
                                            end
                                            # Note, in isolation_array, col 1 corresponds to time 0,
                                                #                           col 2 corresponds to time 1 etc
                                                # So array_time_idx = time + 1, populates array
                                                # in column corresponding to day value of "time"
                                        end

                                        # Check current isolation time does not already exceed the specified
                                        # duration of symptomatic isolation
                                        # If it does, then individual does not enter isolation
                                        if time_left_in_symp_isol > 0
                                            # Supercedes household isolation & contact traced isolation
                                            # Reset that status to zero for all remaining days
                                            for time_idx = time:endtime
                                                array_time_idx = time_idx + 1
                                                states.hh_in_isolation_array[node_itr,array_time_idx] = 0
                                                states.CT_isolation_array[node_itr,array_time_idx] = 0
                                            end
                                        end
                                    end

                                    # Get time since index case first displayed symptoms
                                    index_case_timesymp = states.timesymp[node_itr]
                                    time_left_in_hh_isol = household_isoltime - index_case_timesymp
                                        # Subtract 1, as states.timesymp[node_itr] = 1 corresponds to "0" days elapsed etc

                                    # Irrespective of whether index case self-isolates,
                                    # adherent members of their household may also isolate.
                                    if time_left_in_hh_isol > 0
                                        for hh = 1:household_contacts_per_node[node_itr]
                                            contact_ID = household_contacts[node_itr][hh]
                                            if ( (states.hh_isolation[contact_ID]==1) &&
                                                 (states.symp_isolation_array[contact_ID,output_time_idx]==0) ) # Household member not already symptomatic themselves

                                                # Household member enters household isolation
                                                # Populate household isolation time tracker array
                                                for time_itr = 1:time_left_in_hh_isol
                                                    array_time_idx = time + time_itr
                                                    if array_time_idx <= (endtime + 1)
                                                       states.hh_in_isolation_array[contact_ID,array_time_idx] = 1
                                                    end
                                                    # Note, in isolation_array, col 1 corresponds to time 0,
                                                        #                           col 2 corresponds to time 1 etc
                                                        # So array_time_idx = time + 1, populates array
                                                        # in column corresponding to day value of "time"
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                """
                Assign outputs
                """
                # Assign counts in each disease state to array
                output.numlat[output_time_idx,count,intervention_set_itr] = output.numlat[output_time_idx-1,count,intervention_set_itr]
                output.numinf[output_time_idx,count,intervention_set_itr] = output.numinf[output_time_idx-1,count,intervention_set_itr]
                output.numrep[output_time_idx,count,intervention_set_itr] = output.numrep[output_time_idx-1,count,intervention_set_itr]

                """
                Increment counters
                """
                # Increment counters if node is currently in that state.
                increment_counters!(states) # in additional_fns.jl

                """
                Increment infection process
                """
               # If come to the end of latent time, move to infectious time etc
               for node_itr = 1:cmax
                    # if the node has reached the end of latent infection
                    if states.timelat[node_itr]>states.lattime[node_itr]
                        # move to being infectious
                        states.timelat[node_itr] = -1
                        states.timeinf[node_itr] = 1

                        # Increment time series counts
                        output.numinf[output_time_idx,count,intervention_set_itr] += 1
                        output.newinf[output_time_idx,count,intervention_set_itr] += 1

                        # check if new infected will be asymptomatic
                        if states.asymp[node_itr] > 0
                            output.newasymp[output_time_idx,count,intervention_set_itr] += 1
                        end

                        # check if it is worker that is newly infected
                        if worker_nodes[node_itr].returned_to_work==1
                            output.workersinf[output_time_idx,count,intervention_set_itr] += 1

                            # Count workers that are asymptomatically infected
                            if states.asymp[node_itr] > 0
                                output.workersasymp[output_time_idx,count,intervention_set_itr] += 1
                            elseif workplace_closure_active==true
                                # if the newly infected worker is symptomatic, add to
                                # the workplace memory
                                workertype_ID = worker_nodes[node_itr].sector_ID
                                workplace_ID = worker_nodes[node_itr].workplace_ID
                                workplace_memory[workertype_ID][workplace_ID,WP_memory_slot] += 1
                            end
                        end
                    end

                    # Update node disease state time vectors
                    if states.timeinf[node_itr]>states.inftime
                        # the node becomes symptomatic (if they develop symptoms)
                        states.timeinf[node_itr] = -1
                        states.timesymp[node_itr] = 1

                        # Increment time series counts
                        output.numrep[output_time_idx,count,intervention_set_itr] += 1

                        # Check if index case are symptomatic & would have zero adherence delay
                        if (states.asymp[node_itr] == 0) && (states.delay_adherence[node_itr]==0)
                            # Check if infected will isolate
                            if (states.hh_isolation[node_itr]==1)

                                # Set that the unit has reported infection this timestep
                                states.rep_inf_this_timestep[node_itr] = 1

                                # If testing active,
                                # check if test result will be a false negative
                                if contact_tracing_active == true
                                    tot_time_inf = time - states.acquired_infection[node_itr]

                                    # Get relevant test sensitivity value based on time since infection
                                    if tot_time_inf == 0
                                        test_detection_prob = 0.
                                    else
                                        test_detection_prob = test_detection_prob_vec[tot_time_inf]
                                    end

                                    # Bernoulli trial to determine if false negative returned
                                    if r_symp_test[node_itr,time] < (1 - test_detection_prob)
                                        CT_vars.Test_result_false_negative[node_itr] = true
                                    end
                                end

                                # If testing taking place, check if test result is positive or negative
                                # Set length of isolation
                                if CT_vars.Test_result_false_negative[node_itr] == true
                                   # Release from symptomatic isolation once test result received
                                   # (if does not exceed usual symptomatic isolation period)
                                   length_of_symp_isol = min(CT_vars.CT_delay_until_test_result[node_itr],symp_isoltime)
                                else
                                   # Full symptomatic period spent in isolation
                                   length_of_symp_isol = symp_isoltime

                                   # Result will return a positive.
                                   # Supercedes household isolation & contact traced isolation
                                   # Reset that status to zero for all remaining days
                                   for time_idx = time:endtime
                                       array_time_idx = time_idx + 1
                                       states.hh_in_isolation_array[node_itr,array_time_idx] = 0
                                       states.CT_isolation_array[node_itr,array_time_idx] = 0
                                   end
                                end

                                # Isolation due to presence of symptoms
                                for time_itr = 1:length_of_symp_isol
                                    array_time_idx = time + time_itr
                                    if array_time_idx <= (endtime + 1)
                                        states.symp_isolation_array[node_itr,array_time_idx] = 1
                                    end
                                        # Note, in isolation_array, col 1 corresponds to time 0,
                                        #                           col 2 corresponds to time 1 etc
                                        # So array_time_idx = time + 1, populates array
                                        # in column corresponding to day value of "time"
                                end
                            end

                            # If testing taking place, check if test result is positive or negative
                            # Set length of isolation
                            if CT_vars.Test_result_false_negative[node_itr] == true
                               # Test result will be a negative. Set isolation length
                               # to be used by adherent hh members to result delay time
                               length_of_hh_isol = min(CT_vars.CT_delay_until_test_result[node_itr],household_isoltime)
                            else
                               # Household members to spend full time in isolation
                               length_of_hh_isol = household_isoltime
                            end

                            # Irrespective of whether index case self-isolates,
                            # adherent members of their household may also isolate.
                            for hh = 1:household_contacts_per_node[node_itr]
                                contact_ID = household_contacts[node_itr][hh]
                                if (states.hh_isolation[contact_ID]==1) &&
                                    (states.symp_isolation_array[contact_ID,output_time_idx]==0)  # Household member not already symptomatic themselves

                                    # Populate household isolation time tracker array
                                    for time_itr = 1:length_of_hh_isol
                                        array_time_idx = time + time_itr
                                        if array_time_idx <= endtime+1
                                                states.hh_in_isolation_array[contact_ID,array_time_idx] = 1
                                        end
                                            # Note, in isolation_array, col 1 corresponds to time 0,
                                            #                           col 2 corresponds to time 1 etc
                                            # So array_time_idx = time + 1, populates array
                                            # in column corresponding to day value of "time"
                                    end
                                end
                            end
                        end
                    end

                    # Check if node, if having a delayed adherence, begins adherence on current day
                    if (states.timesymp[node_itr] > 1)&&((states.timesymp[node_itr]-1)==states.delay_adherence[node_itr]) # Condition for node beginning adherence on current day & has been symptomatic for at least one day
                        if states.asymp[node_itr] == 0 # Check node is symptomatic and will adhere
                            if states.hh_isolation[node_itr]==1 # Check node will adhere

                                # Set that the unit has reported infection this timestep
                                states.rep_inf_this_timestep[node_itr] = 1

                                # If testing active,
                                # check if test result will be a false negative
                                if contact_tracing_active == true
                                    tot_time_inf = time - states.acquired_infection[node_itr]

                                    # Get relevant test sensitivity value based on time since infection
                                    if tot_time_inf == 0
                                        test_detection_prob = 0.
                                    else
                                        test_detection_prob = test_detection_prob_vec[tot_time_inf]
                                    end

                                    # Bernoulli trial to determine if false negative returned
                                    if r_symp_test[node_itr,time] < (1 - test_detection_prob)
                                        CT_vars.Test_result_false_negative[node_itr] = true
                                    end
                                end

                                # If testing taking place, check if test result is positive or negative
                                # Set length of isolation.
                                # Individual shortens isolation by length of time since
                                # unwell individual began displaying symptoms.
                                if CT_vars.Test_result_false_negative[node_itr] == true
                                    # Release from symptomatic isolation once test result received
                                    # (if does not exceed remainder of symptomatic isolation period)
                                    length_of_symp_isol = min(CT_vars.CT_delay_until_test_result[node_itr],
                                                                symp_isoltime - states.delay_adherence[node_itr])
                                else
                                   # Rest of symptomatic period spent in isolation
                                   length_of_symp_isol = max(0,symp_isoltime - states.delay_adherence[node_itr])

                                   # Result will return a positive.
                                   # Supercedes household isolation & contact traced isolation
                                   # Reset that status to zero for all remaining days
                                   for time_idx = time:endtime
                                       array_time_idx = time_idx + 1
                                       states.hh_in_isolation_array[node_itr,array_time_idx] = 0
                                       states.CT_isolation_array[node_itr,array_time_idx] = 0
                                   end
                                end

                                # Isolation due to presence of symptoms
                                for time_itr = 1:length_of_symp_isol
                                        # Shorten the isolation period by time already elapsed since
                                        # symptom onset
                                    array_time_idx = time + time_itr
                                    if array_time_idx <= (endtime + 1)
                                        states.symp_isolation_array[node_itr,array_time_idx] = 1
                                    end
                                end
                            end

                            # If testing taking place, check if test result is positive or negative
                            # Set length of isolation.
                            # Individual shortens isolation by length of time since
                            # unwell individual began displaying symptoms.
                            if CT_vars.Test_result_false_negative[node_itr] == true
                               # Test result will be a negative. Set isolation length
                               # to be used by adherent hh members to result delay time
                               length_of_hh_isol = min(CT_vars.CT_delay_until_test_result[node_itr],
                                                    household_isoltime - states.delay_adherence[node_itr])
                            else
                               # Household members to spend full time in isolation
                               length_of_hh_isol = max(0,household_isoltime - states.delay_adherence[node_itr])
                            end

                            # Irrespective of whether index case self-isolates,
                            # adherent members of their household may also isolate.
                            for hh = 1:household_contacts_per_node[node_itr]
                                contact_ID = household_contacts[node_itr][hh]
                                if (states.hh_isolation[contact_ID]==1) &&
                                    (states.symp_isolation_array[contact_ID,output_time_idx]==0) # Household member not already symptomatic themselves

                                    # Populate household isolation time tracker array
                                    for time_itr = 1:length_of_hh_isol
                                        array_time_idx = time + time_itr
                                        if array_time_idx <= (endtime + 1)
                                            states.hh_in_isolation_array[contact_ID,array_time_idx] = 1
                                        end
                                            # Note, in isolation_array, col 1 corresponds to time 0,
                                            #                           col 2 corresponds to time 1 etc
                                            # So array_time_idx = time + 1, populates array
                                            # in column corresponding to day value of "time"
                                    end
                                end
                            end
                        end
                    end

                    # Check if node has reached end of symptom period
                    if states.timesymp[node_itr]>states.symptime
                        states.timesymp[node_itr] = -1
                    end
                end

                """
                Get daily isolation & atwork status of nodes
                """
                # record whether nodes are in isolation
                for node_itr = 1:cmax
                    # Now record whether node is in isolation on current day
                    if (isolation>0) &&
                        ((states.hh_in_isolation_array[node_itr,output_time_idx] == 1) ||
                        (states.symp_isolation_array[node_itr,output_time_idx] == 1) ||
                        (states.CT_isolation_array[node_itr,output_time_idx] == 1))
                        # Default value is 0. So only update if node is in isolation for any reason
                        contacts.daily_record_inisol[time,node_itr] = 1
                    end

                    # Record whether node is at workplace on current day
                    # Needs to be returned to work, and a day where at workplace
                    # AND not in isolation
                    # AND whole workplace not closed
                    workertype_ID = worker_nodes[node_itr].sector_ID
                    workplace_ID = worker_nodes[node_itr].workplace_ID
                    node_workplace_info = workplace_info[workertype_ID][workplace_ID]
                    if (worker_nodes[node_itr].returned_to_work==1) &&
                        (states.atwork[node_itr,time] == true) &&
                        (contacts.daily_record_inisol[time,node_itr] == false) &&
                        (node_workplace_info.workplace_open == true)

                        # Default value is 0. So only update if node is at workplace
                        contacts.daily_record_atworkplace[time,node_itr] = 1
                    end
                end

                """
                Transmit infections

                Structure:
                - Household
                - At work
                    -- Social contacts
                    -- Work contacts (with checks based on covid-secure status)
                    -- Dynamic contacts
                - Not at work
                """
                # Iterate over nodes that may be able to transmit infection
                for node_itr = 1:cmax
                    if ((states.timeinf[node_itr]>0) | (states.timesymp[node_itr]>0))
                            # Only enter loop if node is capable of transmitting infection

                        # find the total time infectious
                        if states.timeinf[node_itr]>0
                            tot_time_infectious = states.timeinf[node_itr]
                        else
                            tot_time_infectious = states.timesymp[node_itr]+states.inftime
                        end

                        # find the infectiousness
                        infectiousness = dist_infectivity[tot_time_infectious]
                        current_worker = worker_nodes[node_itr]
                        if states.asymp[node_itr]>0 # Asymptomatic
                            transtemp_household = current_worker.transrisk_household*infectiousness*asymp_trans_scaling
                            transtemp_work_static = current_worker.transrisk_static_work*infectiousness*asymp_trans_scaling
                            transtemp_work_dynamic = current_worker.transrisk_dynamic_work*infectiousness*asymp_trans_scaling
                            transtemp_social = current_worker.transrisk_social*infectiousness*asymp_trans_scaling
                            transtemp_random = current_worker.transrisk_random*infectiousness*asymp_trans_scaling
                        else
                            if states.timesymp[node_itr]>0  # symptomatic & less infectious due to cautionary behaviour
                                transtemp_household = current_worker.transrisk_household*infectiousness*iso_trans_scaling
                                transtemp_work_static = current_worker.transrisk_static_work*infectiousness*iso_trans_scaling
                                transtemp_work_dynamic = current_worker.transrisk_dynamic_work*infectiousness*iso_trans_scaling
                                transtemp_social = current_worker.transrisk_social*infectiousness*iso_trans_scaling
                                transtemp_random = current_worker.transrisk_random*infectiousness*iso_trans_scaling
                            else  # infected, unscaled infectiousness
                                transtemp_household = current_worker.transrisk_household*infectiousness
                                transtemp_work_static = current_worker.transrisk_static_work*infectiousness
                                transtemp_work_dynamic = current_worker.transrisk_dynamic_work*infectiousness
                                transtemp_social = current_worker.transrisk_social*infectiousness
                                transtemp_random = current_worker.transrisk_random*infectiousness
                            end
                        end

                        # Infection check for other household members
                        # Transmit over household_contacts[node_itr]
                        # checking that contacts are susceptible
                        n_hh_contacts = length(household_contacts[node_itr])
                        if n_hh_contacts > 0
                               transmit_over!(transtemp_household,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                    infecting_by = node_itr,
                                                    contacts_to_check = household_contacts[node_itr],
                                                    transmission_setting="household")
                        end

                        # Check individual is not in household isolation.
                        # If not, can see if transmission across cohort, society,
                        # dynamic household contacts occured.
                        if (contacts.daily_record_inisol[time,node_itr] == false)
                            # Check if node is at workplace or not
                            if contacts.daily_record_atworkplace[time,node_itr] == true
                                if contacts.social_contacts_per_node[node_itr] > 0
                                    # Satisfied condition that node_itr may have social links
                                    # transmit over contacts.workday_social_contacts_by_day[time,node_itr]
                                    # checking that contacts are susceptible and not isolating
                                    transmit_over!(transtemp_social,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                            infecting_by=node_itr,
                                            contacts_to_check=contacts.workday_social_contacts_by_day[time,node_itr],
                                            inisol=contacts.daily_record_inisol,
                                            atwork=undefined_array,
                                            transmission_setting="social",
                                            social_contact_scaling=network_parameters.social_contact_scaling[node_itr])
                                end

                                # Check if workplace is Covid-Secure
                                # If so, modify the transmission risk and use CS contacts
                                # Also, non-workplace worker contacts made
                                workertype_ID = current_worker.sector_ID
                                workplace_ID = current_worker.workplace_ID
                                current_workplace_info = network_parameters.workplace_info[workertype_ID][workplace_ID]
                                if (CS_active_flag == true) && (current_workplace_info.covid_secure == true)
                                    transtemp_work_static *= CS_scale_transrisk[workertype_ID]
                                    transtemp_work_dynamic *= CS_scale_transrisk[workertype_ID]
                                    # These are the regular worker contacts in the same workplace
                                    # in a CS setting
                                    n_work_contacts_same_workplace_CS = contacts.work_contacts_same_workplace_per_node_CS[node_itr]
                                    if n_work_contacts_same_workplace_CS > 0
                                        # Satisfied condition that node_itr has any workday links
                                        # transmit over contacts.work_contacts_same_workplace_CS[time,node_itr]
                                        # checking that contacts are susceptible, not isolating and at work
                                        transmit_over!(transtemp_work_static,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                infecting_by = node_itr,
                                                contacts_to_check = contacts.work_contacts_same_workplace_CS[node_itr],
                                                inisol = contacts.daily_record_inisol,
                                                atwork = contacts.daily_record_atworkplace,
                                                transmission_setting="work")
                                    end
                                else
                                    # These are the regular worker contacts in the same workplace
                                    n_work_contacts_same_workplace = work_contacts_same_workplace_per_node[node_itr]
                                    if n_work_contacts_same_workplace > 0
                                        # Satisfied condition that node_itr has any workday links
                                        # transmit over contacts.work_contacts_same_workplace[node_itr]
                                        # checking that contacts are susceptible, not isolating and at work
                                        transmit_over!(transtemp_work_static,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                infecting_by = node_itr,
                                                contacts_to_check = contacts.work_contacts_same_workplace[node_itr],
                                                inisol = contacts.daily_record_inisol,
                                                atwork = contacts.daily_record_atworkplace,
                                                transmission_setting="work")
                                    end

                                    # These are the regular worker contacts with workers from
                                    # other workplaces
                                    n_work_contacts_other_workplace = work_contacts_other_workplace_per_node[node_itr]
                                    if n_work_contacts_other_workplace > 0
                                        # Satisfied condition that node_itr has any workday links
                                        # transmit over contacts.work_contacts_other_workplace[node_itr]
                                        # checking that contacts are susceptible, not isolating and at work
                                        # also check their workplace is not CS - if CS then no close contact occurs.
                                        transmit_over_other_workplace!(transtemp_work_static,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                                        node_itr,
                                                                        contacts.work_contacts_other_workplace[node_itr],
                                                                        contacts.daily_record_inisol,
                                                                        contacts.daily_record_atworkplace,
                                                                        network_parameters)
                                    end
                                end

                                # Add dynamic links, if returned to work
                                # transmit over contacts.dynamic_worker_contacts[time,node_itr]
                                # checking that contacts are susceptible and not isolating
                                transmit_over!(transtemp_work_dynamic,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                        infecting_by = node_itr,
                                        contacts_to_check = contacts.dynamic_worker_contacts[time,node_itr],
                                        inisol = contacts.daily_record_inisol,
                                        atwork=undefined_array,
                                        dynamic_contact = 1,
                                        transmission_setting="work")

                            else # otherwise the node is not at work

                                # Add in social contacts on non-work days
                                if contacts.social_contacts_per_node[node_itr] > 0
                                    # Satisfied condition that node_itr may have social links
                                    # transmit over contacts.nonworkday_social_contacts_by_day[time,node_itr]
                                    # checking that contacts are susceptible and not isolating
                                    transmit_over!(transtemp_social,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                            infecting_by=node_itr,
                                            contacts_to_check=contacts.nonworkday_social_contacts_by_day[time,node_itr],
                                            inisol=contacts.daily_record_inisol,
                                            atwork=undefined_array,
                                            transmission_setting="social",
                                            social_contact_scaling=network_parameters.social_contact_scaling[node_itr])
                                end
                            end

                            # If not in isolation, also include random daily contacts
                            n_random_contacts = length(dynamic_random_contacts[time, node_itr])
                            if n_random_contacts > 0
                               transmit_over!(transtemp_random,states.timelat,infected_by,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                    infecting_by = node_itr,
                                                    contacts_to_check = dynamic_random_contacts[time, node_itr],
                                                    transmission_setting="other",
                                                    random_contact_scaling=network_parameters.random_contact_scaling[node_itr])
                            end
                        end
                    end
                end

                """
                Perform contact tracing
                """
                # If in use, enact contact tracing from index cases reported today
                if contact_tracing_active == true

                        # Store contacts made during day.
                        # For those reporting symptoms, start delay to test result (if needed)
                        for node_itr = 1:cmax

                            # Increment time to test result if currently waiting for that to occur
                            if CT_vars.Time_to_test_result[node_itr]>=0
                                CT_vars.Time_to_test_result[node_itr] += 1
                            end

                            # For current worker, check if would be leaving pre-symptomatic phase
                            # and displaying symptoms
                            # If so, and will not return a false negative test result, gather traceable contacts
                            if (states.rep_inf_this_timestep[node_itr] == 1)

                                # Increment test counter
                                output.tests_performed[output_time_idx,count,intervention_set_itr] += 1

                                # Initialise CT_vars.Time_to_test_result value
                                CT_vars.Time_to_test_result[node_itr] = 0

                                # Check if the worker is destined to return a false negative result
                                # If false negative to be returned, do not need to work out who traceable contacts are
                                # Otherwise, gather traceable contacts
                                # Also, if no isolation in use, no need to gather contacts.
                                CT_vars.Inds_to_be_contacted[node_itr] = Int64[] # Initialise vector to store contacts
                                if (CT_vars.Test_result_false_negative[node_itr] == false) &&
                                        (isolation>0) && (CT_vars.Engage_with_CT[node_itr] == true)

                                    trace_node!(node_itr,time,CT_vars,contacts,CT_parameters,network_parameters,rng)

                                    # if we are doing "backward contact tracing"
                                    # some small chance the infector is included in this
                                    # don't try to backwards trace the initial infections
                                    if (r_backwards_CT[node_itr,time]<prob_backwards_CT) && (infected_by[node_itr]!=-1)
                                        append!(CT_vars.Inds_to_be_contacted[node_itr],infected_by[node_itr])
                                            CT_vars.Recall_infector[node_itr] = 1
                                    end
                                end

                                # Remove duplicates in CT_vars.Inds_to_be_contacted[node_itr]
                                unique!(CT_vars.Inds_to_be_contacted[node_itr])
                            end

                            # Check if delay to test result reached
                            if CT_vars.Time_to_test_result[node_itr] >= CT_vars.CT_delay_until_test_result[node_itr]

                                # Reset the CT_vars.Time_to_test_result counter
                                CT_vars.Time_to_test_result[node_itr] = -1

                                # If delay time passed, check if test returned a false negative.
                                if CT_vars.Test_result_false_negative[node_itr] == true

                                    # Increment false negative counter
                                    output.test_outcomes[output_time_idx,count,intervention_set_itr,2] += 1

                                    # # Amend tracker of symptomatic cases, unknown test result
                                    # current_node_household_ID = worker_nodes[node_itr].household_ID
                                    # CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1

                                    # # Error check
                                    # if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
                                    #     error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
                                    # end
                                else
                                    # Increment true positive counter
                                    output.test_outcomes[output_time_idx,count,intervention_set_itr,1] += 1

                                    # If test is positive, contacts told to isolate
                                    # Get number of recallable contacts
                                    n_recallable_contacts = length(CT_vars.Inds_to_be_contacted[node_itr])
                                    output.num_CT[output_time_idx,count,intervention_set_itr] += n_recallable_contacts

                                    # Do contacts isolate or not?
                                    # Isolation based on adherence to household isolation of that individual
                                    if isolation>0
                                        for recallable_contact_idx = 1:n_recallable_contacts
                                            recallable_contact_ID = CT_vars.Inds_to_be_contacted[node_itr][recallable_contact_idx]

                                            # Check if individual will adhere to guidance
                                            # If so, they self-isolate
                                            if states.hh_isolation[recallable_contact_ID] == 1
                                                # Populate contact tracing isolation time tracker array
                                                # Time needed to spend in isolation reduced by test result delay time for index case
                                                time_in_CT_isol = max(0,CT_caused_isol_limit - CT_vars.CT_delay_until_test_result[node_itr])
                                                for time_itr = 1:time_in_CT_isol
                                                    array_time_idx = time + time_itr
                                                    if array_time_idx <= (endtime + 1)
                                                        states.CT_isolation_array[recallable_contact_ID,array_time_idx] = 1
                                                    end
                                                        # Note, in isolation_array, col 1 corresponds to time 0,
                                                        #                           col 2 corresponds to time 1 etc
                                                        # So array_time_idx = time + 1, populates array
                                                        # in column corresponding to day value of "time"
                                                end
                                            end
                                        end

                                        # Perform forwards CT from infector, if infector has been identified
                                        # and such a policy is active
                                        if perform_CT_from_infector == true
                                            if CT_vars.Recall_infector[node_itr] == 1
                                                recall_infector_count += 1

                                                forwardCT_from_infector!(infected_by[node_itr],
                                                                            CT_vars,
                                                                            contacts,
                                                                            CT_parameters,
                                                                            CT_vars.Engage_with_CT,
                                                                            states.atwork,
                                                                            time,
                                                                            count,
                                                                            states.hh_isolation,                                                                            network_parameters,
                                                                            rng,
                                                                            infector_trace_count)
                                            end
                                        end
                                    end
                                end
                            end
                        end
                end

                """
                Assign prevalence & isolation outputs
                """
                # For this timestep, get number isolating
                # and if they are latent infected or infectious on current timestep
                for node_itr=1:cmax
                    # Initialise isolation status flag
                    isolating_for_any_reason = false

                    # Isolating due to housemate symptoms
                    if states.hh_in_isolation_array[node_itr,output_time_idx] == 1
                        output.num_household_isolating[output_time_idx,count,intervention_set_itr] += 1
                        isolating_for_any_reason = true
                    end

                    # Isolating due to symptoms
                    if states.symp_isolation_array[node_itr,output_time_idx] == 1
                        output.num_symp_isolating[output_time_idx,count,intervention_set_itr] += 1
                        isolating_for_any_reason = true
                    end

                    # Isolating as close contact of positive test symptomtic
                    if contact_tracing_active == true
                        if states.CT_isolation_array[node_itr,output_time_idx] == 1
                            output.num_isolating_CTcause[output_time_idx,count,intervention_set_itr] += 1
                            isolating_for_any_reason = true
                        end
                    end

                    # Isolating for any reason
                    if isolating_for_any_reason == true
                        output.num_isolating[output_time_idx,count,intervention_set_itr] += 1
                    end

                    # Check if latently infected
                    if states.timelat[node_itr]>0
                        output.prevlat[output_time_idx,count,intervention_set_itr] += 1
                    end

                    # In presymptomatic infectious period.
                    if (states.timeinf[node_itr]>0)
                        if states.asymp[node_itr] > 0 # asymptomatic
                            output.prevasymp[output_time_idx,count,intervention_set_itr] += 1
                        else # will be symptomatic
                            output.prevpresymp[output_time_idx,count,intervention_set_itr] += 1
                        end
                    end

                    # After presymp period, check if symptomatic or asymptomatic
                    if (states.timesymp[node_itr]>0)
                        if states.asymp[node_itr] > 0 # asymptomatic
                            output.prevasymp[output_time_idx,count,intervention_set_itr] += 1
                        else # symptomatic
                            output.prevsymp[output_time_idx,count,intervention_set_itr] += 1
                        end
                    end

                    # Check if recovered
                    if states.timesymp[node_itr] == -1
                        output.prevrec[output_time_idx,count,intervention_set_itr] += 1
                    end
                end

                """
                Reactive workplace closure check
                """
                # close workplaces with too many infections
                if workplace_closure_active==true
                    for worktypeID = 1:workertypes
                        for work_ID = 1:num_workplaces[worktypeID]
                            if work_closed_time[worktypeID][work_ID]>0
                                # if the workplace is closed, move on the time counter
                                work_closed_time[worktypeID][work_ID] += 1
                            else
                                # otherwise decide if the workplace should be closed
                                total_infections = 0
                                for ii=1:workplace_CT_memory
                                    total_infections += workplace_memory[worktypeID][work_ID,ii]
                                end

                                # Number of infections in time window exceeds threshold
                                # Set workplace to be closed
                                if total_infections>workplace_thresholds[worktypeID][work_ID]
                                    work_closed_time[worktypeID][work_ID] = 1
                                    network_parameters.workplace_info[worktypeID][work_ID].workplace_open = false
                                end
                            end

                            if work_closed_time[worktypeID][work_ID]>time_WC
                                # if the workplace has been closed long enough, open it
                                work_closed_time[worktypeID][work_ID] = 0

                                # Check sector is not closed. If sector is open, workplace set to be open
                                if (network_parameters.sector_open[worktypeID] == true)
                                    network_parameters.workplace_info[worktypeID][work_ID].workplace_open = true
                                end
                            end
                        end
                    end

                end

                """
                Run interventions
                """
                # Check if any interventions are triggered
                # Update statuses as needed
                if isassigned(intervention_fns) # Check if any intervetion were specified

                    # Package health outcome measures that may be used in decision
                    # to enact an intervention
                    intervention_trigger_input_data = intervention_data_feeds(rep_inf_this_timestep = states.rep_inf_this_timestep,
                                                                                numlat = output.numlat[output_time_idx,count,intervention_set_itr],
                                                                                numinf = output.numinf[output_time_idx,count,intervention_set_itr],
                                                                                numrep = output.numrep[output_time_idx,count,intervention_set_itr],
                                                                                newinf = output.newinf[output_time_idx,count,intervention_set_itr])

                    for interv_trig_itr = 1:n_intervention_fns
                        # Check if condition is met
                        condition_fn = intervention_fns[interv_trig_itr,1]
                        eval_condition_fn = condition_fn(intervention_trigger_input_data,
                                                          time,
                                                          network_parameters)
                        if eval_condition_fn == true
                            # If condition satisfied, apply the affect
                            chosen_affect_fn = intervention_fns[interv_trig_itr,2]
                            chosen_affect_fn(network_parameters)
                        end
                    end
                end
            end

            # Find how many nodes were infected by the initial nodes
            initial_nodes = findall(infected_by.==-1)
            sum_infections = 0
            output.num_init_infected[count] = zeros(Int64,length(initial_nodes),n_intervention_sets) # Initialise output array
            for initial_node_it = 1:length(initial_nodes)
                output.num_init_infected[count][initial_node_it,intervention_set_itr] = output.num_infected[initial_nodes[initial_node_it],count,intervention_set_itr]
                sum_infections+=output.num_init_infected[count][initial_node_it,intervention_set_itr]
            end

            # find mean generation time
            if sum_infections>0
                output.mean_init_generation_time[count,intervention_set_itr] = output.mean_init_generation_time[count,intervention_set_itr]/sum_infections
            end

            # divide number of infections by number of infectors to find Rt
            for time=1:(endtime+1)
                # divide by the number of nodes that were infected (entered the latent state)
                # at time
                if time == 1
                    output.Rt[time,count,intervention_set_itr] = output.Rt[time,count,intervention_set_itr] / output.numlat[time,count,intervention_set_itr]
                else
                    output.Rt[time,count,intervention_set_itr] = output.Rt[time,count,intervention_set_itr] / (output.numlat[time,count,intervention_set_itr]-output.numlat[time-1,count,intervention_set_itr])
                end
            end

            # Print to screen info on run just completed
            println("Run $count complete.")
        end

    end

    # Compute variance in number of infected per node
    # var_num_infected = zeros(countfinal)
    # for count=1:countfinal
    #     var_num_infected[count] = var(output.num_infected[:,count,intervention_set_itr])
    # end

    # Specify what is output from the function
    return output
end
