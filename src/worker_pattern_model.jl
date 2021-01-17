"""
Purpose:
Toy network model to explore impact of working patterns on transmission in the workplace
"""

"""
Set paths & load environment
"""

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../")

"""
Load packages
"""
#Required packages
using MAT, Distributions, CSV
using LinearAlgebra, Random, DelimitedFiles, Parameters

# using Profile  # For use with memory allocation checks

"""
Load supporting files
"""
# Files containing other functions needed to run the model
include("include_files_network_model/parametertypes.jl")
include("include_files_network_model/contact_tracing_fns.jl")
include("include_files_network_model/network_generation_fns.jl")
include("include_files_network_model/workplace_size_generation_fns.jl")
include("include_files_network_model/additional_fns.jl")
include("include_files_network_model/intervention_condition_affect_fns.jl")
include("include_files_network_model/seed_initial_states_fns.jl")
include("include_files_network_model/main_function.jl")

"""
Set variables from ARGS
"""
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] job_ID
# args[2] RNGseed: To be used to initialise the random number generator
# args[3] countfinal: Number of replicates requested
# args[4] endtime: Timesteps for each individual simulation
# args[5] n_nodes: Overall size of network
# args[6] n_sectors: Defining total number of work sectors
# args[7] seed_initial_states_fn
# args[8] runsets: Scenarios to run
if length(ARGS)==0
    args = [ "1", "100", "100", "365", "5000", "41", "seed_states_with_uncertainty", "[\"adherence_intervention\"]"]
end

# To run from command line, example:
# julia worker_pattern_model.jl 1 100 10 365 10000 41 set_ten_initial_infected '["run_one_run"]'

# Example runset options (see load_configs in "include_files_network_model/additional_fns.jl"):
#["synchronised_changedays", "variable_changedays",
#"groups_off", "weeks_on", "amount_backwards_CT",
#"workplace_CT_threshold", "CT_engagement",
#"asymp_scen", "change_scaling",
#"RNGseed_svty","initinfected_svty",
#"transscaling_svty", "adherence_svty",
#"clustering_svty","workplace_clustering_svty",
#"social_clustering_svty","popsize_svty",
#"ER_no_control","CS_workplace_no_control",
# "workplace_clustering_svty_no_social",
# "workplace_clustering_svty_no_household",
# "workplace_clustering_svty_workplace_only",
# "workplace_clustering_svty_no_dynamic",
# "workplace_clustering_svty_workplace_static_only",
# "workplace_clustering_svty_workplace_static_and_social",
# "workplace_clustering_svty_workplace_static_and_household",
# "work_percent_svty",
# "synchronised_changedays_intervention","variable_changedays_intervention",
# "weeks_on_intervention","workpercent_intervention",
# "amount_backwards_CT_intervention", "adherence_intervention"
# "CS_intervention","CS_intervention_no_isol"]

# Set identifier for job
job_ID = parse(Int64, args[1])

# Set RNG seed
RNGseed_base = parse(Int64, args[2])
RNGseed = job_ID*RNGseed_base

# Set simulation run params
countfinal = parse(Int64, args[3])  # Number of simulations to be performed per scenario
endtime = parse(Int64, args[4]) # Timesteps for each individual simulation
cmax= parse(Int64, args[5]) # Number of individuals in the network
workertypes = parse(Int64, args[6]) # Defining total number of work sectors

# # Specify function for setting proportion/number of nodes in each initial disease state
# # Relevant files in "include_files_network_model/seed_initial_states_fn.jl"
# # e.g. set_ten_initial_infected
s = Symbol(args[7])
seed_initial_states_fn = getfield(Main, s) #Make Symbol a callable function

# Specify scenario to be run
runsets = eval(Meta.parse(args[8]))

# Set proportion of nodes to be in recovered state
recov_propn = 0.

"""
Set up incubation & infectivity distributions
"""

# If needed, set up a different latent period distribution
# Default: Erlang(6,0.88)
# d_incub_alt = LogNormal(log(5.2), log(1.4))
    # If using alternative then need to pass d_incub_alt into infection_params creation:
    #e.g. infection_parameters = infection_params(transrisk=0.75*ones(workertypes),
    #                                                  d_incub = d_incub_alt)

# If needed, set a different distribution of infectivity from the default
# Will then need to pass dist_infectivity into infection_params creation:
#  e.g. infection_parameters = infection_params(transrisk=0.75*ones(workertypes),
#                                                  dist_infectivity = dist_infectivity)
# dist_infectivity = ones(10) # if we don't want to have a distribution of infectiousness


"""
Set up adherence delay distribution
"""
# Set a different distribution for delay in adhering in guidance
# Will need to pass delay_adherence_pmf_alt into infection_params creation:
#  e.g. infection_parameters = infection_params(...,
#                                                  delay_adherence_pmf = delay_adherence_pmf_alt)
# delay_adherence_pmf_alt = [0.5,0.25,0.25]
#
#  # Check delay_adherence_pmf_alt sums to 1
# if sum(delay_adherence_pmf_alt) != 1
#     error("delay_adherence_pmf_alt must sum to 1. Currently sums to $(sum(delay_adherence_pmf_alt))")
# end

"""
Set up testing related distributions
"""

# If needed, set a different probability mass function for delay until test result received.
# Will then need to pass CT_delay_until_test_result_pmf_alt into CT_params creation:
#  e.g. CT_parameters = CT_params(...,
#                                  CT_delay_until_test_result_pmf = CT_delay_until_test_result_pmf_alt)
# CT_delay_until_test_result_pmf_alt = [0.,0.,0.5,0.,0.,0.,0.5]
#
#  # Check CT_delay_until_test_result_pmf_alt sums to 1
# if sum(CT_delay_until_test_result_pmf_alt) != 1
#     error("CT_delay_until_test_result_pmf_alt must sum to 1. Currently sums to $(sum(CT_delay_until_test_result_pmf_alt))")
# end

# If needed, set an alternative false negative test prob. , wrt days since infected
# Will need to pass test_false_negative_vec_alt into CT_params creation:
#  e.g. CT_parameters = CT_params(...,
#                                 test_false_negative_vec = test_false_negative_vec_alt)
# test_false_negative_vec_alt = 0.2*ones(20)

"""
Specify use of any additional, trigger interventions
"""
# Have as a 2D array input.
# -> Row per intervention
# -> Column 1 for the condition for itervention being triggered
# -> Column 2 for the affect if condition is satisfied

# intervention_fns_alt = [condition_close_example affect_close_example!;
#                         condition_open_example affect_open_example!]

    # Will need to pass as optional input to worker_pattern_network_run:
    #  e.g. ...=  worker_pattern_network_run(...,
    #                                 intervention_fns = intervention_fns_alt)

"""
Specify household transmission risk by group
"""
# # If needed, set a different distribution of transrisk_household_group from the default
# # Will then need to pass transrisk_household_group_altinto infection_params creation:
# #  e.g. infection_parameters = infection_params(...,
# #                                                  transrisk_household_group = transrisk_household_group_alt)
# transrisk_household_group_alt = [1.,0.5,0.2]
#
# # Specify function to perform household group allocation
assign_household_transrisk_fn_alt = assign_household_transmit_household_size!
#         # Will need to pass as optional input to worker_pattern_network_run:
#         #  e.g. ...=  worker_pattern_network_run(...,
#         #                                 assign_household_transrisk_fn = assign_household_transrisk_fn_alt)

for run_it = 1:length(runsets)

    # Set if contact tracing is active or not (Bool type variable)
    contact_tracing_active = false

    # Set if workplace closures is active or not (Bool type variable)
    workplace_closure_active = false

    # set if backwards contact tracing is active or not (Bool type variable)
    perform_CT_from_infector = false

    runset = runsets[run_it] #  options: "run_one_run", "synchronised_changedays", "variable_changedays", "groups off", "weeks on"

    sameday_config, ton_config, toff_config, work_perc_config,
    num_config, workplace_CT_threshold_config, prob_backwards_CT_config,
    infector_engage_with_CT_prob_config, CT_engagement_config,
    transiso_config, cmax_config,
    rng_config, trans_scaling_config, adherence_config,
    clustering_config, CS_team_size_config, intervention_list_config,
    dynamic_time_frame_config, group_limit_config, max_contacts_social_config = load_configs(runset,workertypes,cmax,RNGseed)

    n_intervention_sets = length(intervention_list_config)

    # Initialise output arrays. Store counts for each network configuration
    # put countfinal in the 2nd place so that we can concatenate variables later
    numlat_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    numinf_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    numrep_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    prevlat_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    prevsymp_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    prevasymp_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    prevpresymp_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    prevrec_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    newinf_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    workersinf_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    workersasymp_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    newasymp_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    num_isolating_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)      # Number isolating on given day
    num_symp_isolating_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config) # Number isolating due to having symptoms
    num_household_isolating_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config) # Number isolating due to having symptoms
    num_CT_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)             # Number contact traced on given day
    num_infected_save = Array{Array{Int64,3}}(undef,num_config) #zeros(Int64,cmax,countfinal,num_config)
    dynamic_infection_count_save = zeros(Int64,endtime+1,countfinal,n_intervention_sets,num_config)
    var_num_infected_save = zeros(Float64,1,countfinal,n_intervention_sets,num_config) # so that the countfinal is in the 2nd place
    mean_init_generation_time_save = zeros(Float64,1,countfinal,n_intervention_sets,num_config) # so that the countfinal is in the 2nd place
    num_isolating_CTcause_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    Rt_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    transmission_setting_save = zeros(Int64, endtime+1, countfinal, n_intervention_sets, 5, num_config)
    tests_performed_save = zeros(endtime+1,countfinal,n_intervention_sets,num_config)
    test_outcomes_save = zeros(endtime+1,countfinal,n_intervention_sets,4,num_config)
    num_init_infected_save = Array{Array{Array{Int64,2}}}(undef,num_config)

    for it = 1:num_config

        sameday = sameday_config[it]
        ton = ton_config[it]
        toff = toff_config[it]
        workpercent = work_perc_config[it,:]
        workplace_CT_threshold = workplace_CT_threshold_config[it]
        prob_backwards_CT = prob_backwards_CT_config[it]
        infector_engage_with_CT_prob = infector_engage_with_CT_prob_config[it]
        CT_engagement = CT_engagement_config[it]
        # asymp_trans_scaling = transasymp_config[it]
        iso_trans_scaling = transiso_config[it]
        global cmax = cmax_config[it]

        ### Sensitivity config setup ###
        RNGseed = rng_config[it]
        scaling_social = trans_scaling_config[it][1]
        scaling_work_static = trans_scaling_config[it][2]
        scaling_work_dynamic = trans_scaling_config[it][3]
        scaling_household = trans_scaling_config[it][4]
        scaling_random = trans_scaling_config[it][5]
        adherence = adherence_config[it]
        CS_team_size = CS_team_size_config[it]

        # Set up array to store number infected by each node
        num_infected_save[it] = zeros(Int64,cmax,countfinal,n_intervention_sets)
        #####

        if runset=="amount_backwards_CT"
            perform_CT_from_infector = true
            contact_tracing_active = true
        end

        if runset=="workplace_CT_threshold"
            workplace_closure_active = true
        end

        if (runset=="CT_engagement")||(runset=="amount_backwards_CT")
            contact_tracing_active = true
        end

        # Set up network and workplace parameters
        # Uses requested number of workertypes
        # For options, the function resides in "additional_fns.jl"
        network_parameters, workplace_generation_parameters = find_network_parameters(workertypes,workpercent=workpercent)
        #network_parameters.CS_active_flag = true


        CT_parameters = CT_params(prob_backwards_CT = prob_backwards_CT,
                                    perform_CT_from_infector = perform_CT_from_infector,
                                    infector_engage_with_CT_prob = infector_engage_with_CT_prob,
                                    workplace_CT_threshold = workplace_CT_threshold,
                                    CT_engagement = CT_engagement,
                                    contact_tracing_active = contact_tracing_active
                                    #CT_delay_until_test_result_pmf = CT_delay_until_test_result_pmf_alt
                                    #test_false_negative_vec = test_false_negative_vec_alt
                                    )

        infection_parameters = infection_params(iso_trans_scaling = iso_trans_scaling,
                                                recov_propn = recov_propn,
                                                #delay_adherence_pmf = delay_adherence_pmf_alt
                                                #transrisk_household_group = transrisk_household_group_alt
                                                )

        #########################################################################
        #### Adjust infection and network parameters for sensitivity configs ####
        #########################################################################
        if runset=="clustering_svty"
            clustering = clustering_config[it]
            network_parameters.between_workplace_contact_probs = [clustering[1]]
            network_parameters.friend_of_friend_prob = clustering[2]
        elseif runset∈["workplace_clustering_svty","workplace_clustering_svty_no_social",
                    "workplace_clustering_svty_no_household","workplace_clustering_svty_workplace_only",
                    "workplace_clustering_svty_no_dynamic","workplace_clustering_svty_workplace_static_only",
                    "workplace_clustering_svty_workplace_static_and_social","workplace_clustering_svty_workplace_static_and_household"]
            clustering = clustering_config[it]
            network_parameters.between_workplace_contact_probs = [clustering]
        elseif runset=="social_clustering_svty"
            clustering = clustering_config[it]
            network_parameters.friend_of_friend_prob = clustering
        end

        if runset=="ER_no_control"
            network_parameters.network_generation_method = "ER"
            network_parameters.network_generation_method_dynamic_social = "ER"
        end

        infection_parameters.transrisk_social_mean *= scaling_social
        infection_parameters.transrisk_static_work_mean *= scaling_work_static
        infection_parameters.transrisk_dynamic_work_mean *= scaling_work_dynamic
        infection_parameters.transrisk_household_group_mean *= scaling_household
        infection_parameters.transrisk_random_mean *= scaling_random

        infection_parameters.transrisk_social_sd *= scaling_social
        infection_parameters.transrisk_static_work_sd *= scaling_work_static
        infection_parameters.transrisk_dynamic_work_sd *= scaling_work_dynamic
        infection_parameters.transrisk_household_group_sd *= scaling_household
        infection_parameters.transrisk_random_sd *= scaling_random

        if runset∈["adherence_svty","synchronised_changedays_intervention",
                    "variable_changedays_intervention","workpercent_intervention",
                    "amount_backwards_CT_intervention","adherence_intervention",
                    "CS_intervention","CS_intervention_no_isol","run_one_run"]
            infection_parameters.adherence = adherence
        end

        if (runset=="CS_workplace_no_control")
            network_parameters.CS_active_flag = true
            network_parameters.CS_team_size = CS_team_size
            infection_parameters.CS_scale_transrisk = ones(workertypes)
        end

        if runset=="workplace_clustering_svty_no_social"
            infection_parameters.transrisk_social_mean *= 0
            infection_parameters.transrisk_social_sd *= 0
        elseif runset=="workplace_clustering_svty_no_household"
            infection_parameters.transrisk_household_group_mean *= 0
            infection_parameters.transrisk_household_group_sd *= 0
        elseif runset=="workplace_clustering_svty_no_dynamic"
            infection_parameters.transrisk_dynamic_work_mean *= 0
            infection_parameters.transrisk_dynamic_work_sd *= 0
        elseif runset=="workplace_clustering_svty_workplace_only"
            infection_parameters.transrisk_social_mean *= 0
            infection_parameters.transrisk_social_sd *= 0
            infection_parameters.transrisk_household_group_mean *= 0
            infection_parameters.transrisk_household_group_sd *= 0
        elseif runset=="workplace_clustering_svty_workplace_static_only"
            infection_parameters.transrisk_social_mean *= 0
            infection_parameters.transrisk_social_sd *= 0
            infection_parameters.transrisk_household_group_mean *= 0
            infection_parameters.transrisk_household_group_sd *= 0
            infection_parameters.transrisk_dynamic_work_mean *= 0
            infection_parameters.transrisk_dynamic_work_sd *= 0
        elseif runset=="workplace_clustering_svty_workplace_static_and_social"
            infection_parameters.transrisk_household_group_mean *= 0
            infection_parameters.transrisk_household_group_sd *= 0
            infection_parameters.transrisk_dynamic_work_mean *= 0
            infection_parameters.transrisk_dynamic_work_sd *= 0
        elseif runset=="workplace_clustering_svty_workplace_static_and_household"
            infection_parameters.transrisk_social_mean *= 0
            infection_parameters.transrisk_social_sd *= 0
            infection_parameters.transrisk_dynamic_work_mean *= 0
            infection_parameters.transrisk_dynamic_work_sd *= 0
        end

        if runset∈["RNGseed_svty","initinfected_svty", "transscaling_svty",
                                    "clustering_svty","popsize_svty","ER_no_control",
                                    "CS_workplace_no_control","workplace_clustering_svty",
                                    "social_clustering_svty","workplace_clustering_svty_no_social",
                                    "workplace_clustering_svty_no_household","workplace_clustering_svty_workplace_only",
                                    "workplace_clustering_svty_workplace_static_and_social",
                                    "workplace_clustering_svty_workplace_static_and_household"]
            infection_parameters.isolation = 0
        end

        ########################################################################

        @time  output = worker_pattern_network_run(RNGseed,
        cmax,
        ton,toff,
        infection_parameters,
        sameday,
        seed_initial_states_fn,
        countfinal,
        endtime,
        contact_tracing_active,
        CT_parameters,
        network_parameters,
        workplace_generation_parameters,
        workplace_closure_active,
        intervention_list_config,
        runset,
        assign_household_transrisk_fn = assign_household_transrisk_fn_alt
        )

        @unpack numlat, numinf, numrep,
        prevlat, prevsymp, prevasymp, prevpresymp, prevrec,
        newinf, newasymp, workersinf, workersasymp, infected_by,
        num_isolating, num_household_isolating, num_symp_isolating, num_isolating_CTcause,
        num_CT, num_infected, dynamic_infection_count,
        var_num_infected, num_init_infected, mean_init_generation_time,
        Rt, transmission_setting, tests_performed, test_outcomes = output

        # For this iteration's network config, write results to output storage arrays
        numlat_save[:,:,:,it] = numlat
        numinf_save[:,:,:,it] = numinf
        numrep_save[:,:,:,it] = numrep
        prevlat_save[:,:,:,it] = prevlat
        prevsymp_save[:,:,:,it] = prevsymp
        prevasymp_save[:,:,:,it] = prevasymp
        prevpresymp_save[:,:,:,it] = prevpresymp
        prevrec_save[:,:,:,it] = prevrec
        newinf_save[:,:,:,it] = newinf
        workersinf_save[:,:,:,it] = workersinf
        workersasymp_save[:,:,:,it] = workersasymp
        newasymp_save[:,:,:,it] = newasymp
        num_isolating_save[:,:,:,it] = num_isolating
        num_symp_isolating_save[:,:,:,it] = num_symp_isolating
        num_household_isolating_save[:,:,:,it] = num_household_isolating
        num_CT_save[:,:,:,it] = num_CT
        num_infected_save[it][:,:,:] = num_infected
        dynamic_infection_count_save[:,:,:,it] = dynamic_infection_count
        var_num_infected_save[1,:,:,it] = var_num_infected
        num_init_infected_save[it] = num_init_infected
        mean_init_generation_time_save[1,:,:,it] = mean_init_generation_time
        Rt_save[:,:,:,it] = Rt
        transmission_setting_save[:,:,:,:,it] = transmission_setting
        tests_performed_save[:,:,:,it] = tests_performed
        test_outcomes_save[:,:,:,:,it] = test_outcomes

        # if contact_tracing_active==true
        num_isolating_CTcause_save[:,:,:,it] = num_isolating_CTcause
        # end
    end

    # Save outputs to file
    output_file = "../Results/worker_model_output_$(runset)_#$(job_ID).mat"
    file = matopen(output_file, "w")
    write(file,"numlat",numlat_save)
    write(file,"numinf",numinf_save)
    write(file,"numrep",numrep_save)
    write(file,"prevlat",prevlat_save)
    write(file,"prevsymp",prevsymp_save)
    write(file,"prevasymp",prevasymp_save)
    write(file,"prevpresymp",prevpresymp_save)
    write(file,"prevrec",prevrec_save)
    write(file,"newinf",newinf_save)
    write(file,"workersinf",workersinf_save)
    write(file,"workersasymp",workersasymp_save)
    write(file,"newasymp",newasymp_save)
    write(file, "num_isolating", num_isolating_save)
    write(file, "num_symp_isolating", num_symp_isolating_save)
    write(file, "num_household_isolating", num_household_isolating_save)
    write(file, "num_infected", num_infected_save)
    write(file, "dynamic_infection_count", dynamic_infection_count_save)
    write(file, "var_num_infected_save", var_num_infected_save)
    write(file, "num_init_infected_save", num_init_infected_save)
    write(file, "mean_init_generation_time_save", mean_init_generation_time_save)
    write(file, "Rt_save", Rt_save)
    write(file, "transmission_setting", transmission_setting_save)
    write(file, "tests_performed", tests_performed_save)
    write(file, "test_outcomes", test_outcomes_save)
    write(file, "num_CT", num_CT_save)
    write(file, "num_isolating_CTcause", num_isolating_CTcause_save)
    if runset=="amount_backwards_CT"
        write(file, "prob_backwards_CT", prob_backwards_CT_config)
    end
    close(file)
end
