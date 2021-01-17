"""
Purpose:
Parameter types to be used with the network model

- workplace_params: Info on, for example, the work sector, being covid secure, whether workplace is open
- worker_params: Worker parameter type to have information on return to work status, the work sector etc
- CTParams: Contact tracing related
- contacts_struct: Used for contact structures, that are also used outside of contact tracing
- contact_tracing_vars: Vars only needed for contact tracing
- infection_params: Parameters related to the transmission process
- network_params: Quantities to construct the contacts & stores the node properties
- workplace_generation_params: Used for generating a "population" of workplaces
- sim_outputs: Structures to hold data to be collected from simulation run
- intervention_data_feeds: The data streams that may be of use in implementing trigger-based interventions
- node_states: Used for the status of each node
- intervention_struct: Variable values to be used from a specified intervention start time
"""

# Workplace parameter type to have information on
# the work sector, being covid secure, whether workplace is open etc
@with_kw mutable struct workplace_params
   covid_secure::Bool = false    # If true, workplace is covid secure. Can alter transmission risk & contacts
   workplace_open::Bool = true   # For imposing workplace closures, can set flag to false
end

# Worker parameter type to have information on return to work status,
# the work sector, and workplace.
# Can add other fields as required
@with_kw mutable struct worker_params
   returned_to_work::Int64                # Flag variable. If 1, worker returns to workplace as designated by atwork schedule.
   sector_ID::Int64                       # The worktype/sector the worker has been allocated to
   workplace_ID::Int64                    # The ID of the workplace (within that worktype/sector)
   household_ID::Int64 = 0                # Household the individual has been assigned to
   transrisk_household::Float64 = 0.      # Secondary attack rate within household were that individual the index case
   transrisk_static_work::Float64 = 0.      # Secondary attack rate within household were that individual the index case
   transrisk_dynamic_work::Float64 = 0.      # Secondary attack rate within household were that individual the index case
   transrisk_social::Float64 = 0.      # Secondary attack rate within household were that individual the index case
   transrisk_random::Float64 = 0.
end

@with_kw mutable struct CT_params
# used for parameters relating to contact tracing

   # flag indicating whether or not CT is active
   contact_tracing_active = false

   # For those adhering to self-isolating, engagement with contact tracing aspect
   # If 1, as well as self-isolating, all adhering individuals do give contacts
   CT_engagement::Float64 = 1.

   # Set time delay. 0 corresponds to result on day of reporting.
   CT_delay_until_test_result_pmf::Array{Float64,1} = [0., 0., 1.,]

   # Set number of days before symptoms CT will attempt to capture
   CT_days_before_symptom_included::Int64 = 2

   # Propotions of tests returning a false negative outcome
   # Entry per day since infected
   test_detection_prob_vec::Array{Float64,1} = [0.,0.11,0.59,0.785,0.83,0.845,0.84,0.82,0.79,0.76,  # Days 1-10
                                                0.72,0.68,0.64,0.59,0.54,0.485,0.445,0.405,0.37,0.335, # Days 11-20
                                                0.30,0.27,0.24,0.22,0.20,0.18,0.16,0.15,0.14,0.13] # Days 21-30

   # Amount of time spent in isolation if contact traced
   CT_caused_isol_limit::Int64 = 14

   # Set up recall of dynamic contacts probability
   # Proportion of contacts remembered x days ago
   dynamic_contacts_recalled_propn::Array{Float64,1} = [0.5,0.4,0.3,0.2,0.1]
   social_contacts_recalled_propn::Array{Float64,1} = [1.,1.,1.,1.,1.]

   # proportion of people that can identify their infector (set to 0 for no backwards CT)
   prob_backwards_CT::Float64 = 0.

   # Parameters for performing forward contact tracing from identified infectors
   perform_CT_from_infector::Bool = false
   infector_engage_with_CT_prob::Float64 = 1.0  # For those not complying with other isolation guidance,
                                       # probability they do if idenfitied as possible infector
                                       # Note, those that are set to adhere to isolation guidance
                                       # are also assuemd fully compliant with CT from infector measures

   # Flag indicating workplace closure active or not
   workplace_closure_active::Bool = false

   workplace_CT_memory::Int64 = 7 # how many days to look over when looking for clusters

   workplace_CT_threshold::Float64 = 0.5 # what proportion of employees need to be infected to trigger workplace closure

   time_WC::Int64 = 14 # how many days to close workplaces for
end

@with_kw mutable struct contacts_struct
# used for contact structures, that are also used outside of contact tracing

   # sizes to initialise arrays
   cmax::Int64 = 0
   endtime::Int64 = 0

   # Stores work contacts at the same workplace
   work_contacts_same_workplace::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax)

   # Stores work contacts at other workplaces
   work_contacts_other_workplace::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax)

   # Store reduced work contacts at the same workplace when COVID secure
   work_contacts_same_workplace_CS::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax)
   work_contacts_same_workplace_per_node_CS::Array{Int64,1} = zeros(Int64,cmax)

   # Flag arrays to keep a log of whether a node is in isolation/at workplace each timestep
   # Used in each replicate
   daily_record_atworkplace::Array{Int64,2} = zeros(Int64,endtime,cmax)
   daily_record_inisol::Array{Int64,2} = zeros(Int64,endtime,cmax)

   # Total social contacts made by node
   social_contacts_per_node::Array{Int64,1} = zeros(Int64,cmax)

   # Per node, a record of social contacts made on each day
   workday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)
   nonworkday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

   # Per node, a record of dynamic worker contacts made on each day
   dynamic_worker_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

   # Random contacts made each day by each node
   dynamic_random_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

end

@with_kw mutable struct contact_tracing_vars
# only needed for contact tracing

   # sizes to initialise arrays
   cmax::Int64 = 0
   endtime::Int64 = 0
   n_households::Int64 = 0

   # The number of days prior to symptoms that each node remembers
   relevant_prev_days_for_CT::Array{Int64,1} = zeros(Int64,cmax)

   # Vector of vectors for storing IDs of those to be contacted in CT
   Inds_to_be_contacted::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax)

   # Vector tracking the latest isolation release time due to household member cases
   hh_isolation_release_time::Array{Int64,1} = zeros(Int64,cmax)

   # Variables for waiting for test results, set at -1 until activated
   Time_to_test_result::Array{Int64,1} = -1*ones(Int64,cmax)

   # Boolean vector to store whether a false negative test result would be returned
   # Populated in each replicate
   Test_result_false_negative::Array{Bool,1} = Array{Bool,1}(undef,cmax)

   # Boolean vector to store if individual will engage with CT or not
   # Populated in each replicate
   Engage_with_CT::Array{Bool,1} = Array{Bool,1}(undef,cmax)

   # Array to keep track of whether an infected recalls their infector
   Recall_infector::Array{Int64,1} = zeros(Int64,cmax)

   # Delay before test result is returned
   CT_delay_until_test_result::Array{Int64,1} = zeros(Int64,cmax)
end

@with_kw mutable struct infection_params
# used for parameters relating to the infection process
# NOTE THAT MANY OF THESE ARE RESET IN THE CONFIGURATION

   # Transmission scaling for those symptomatic (to capture cautionary behaviour if symptomatic)
   iso_trans_scaling::Float64 = 1.

   # Also add in the possibility that individuals are asymptomatic.
   # In that case, their infection potential is lower but they dont isolate.
   asymp_trans_scaling_dist::Uniform{Float64} = Uniform(0.3,0.7)

   # Household size distribution (2 - 5+ only, for use in household transrisk average)
   household_size_distribution::Array{Float64,1} = [0.423, 0.075, 0.025, 0.008]/sum([0.423, 0.075, 0.025, 0.006, 0.002])

   # probability of transmission within a household, based on secondary attack rate
   # Draw from a normal distribution.
   # Can differ per household group.
   transrisk_household_group_mean::Array{Float64,1} = [0.48,0.40,0.33,0.22]
   transrisk_household_group_sd::Array{Float64,1} = [0.06,0.06,0.05,0.05]

   # Averages of household transmission across all sizes
   transrisk_household_mean_average::Float64 = 0.37 #sum(transrisk_household_group_mean.*household_size_distribution)
   transrisk_household_sd_average::Float64 = 0.03 #sum(transrisk_household_group_sd.*household_size_distribution)

   # Baseline risk in workplace. Sector dependent, maximum value is 1.0, corresponding to
   # a risk of transmission given by scaling
   # Drawn from normal distribution - means and sds are relative to household transmission = 0.5
   transrisk_static_work_mean::Array{Float64,1} = (transrisk_household_mean_average/0.5).*[0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.21, 0.21, 0.21, 0.30, 0.30, 0.17, 0.27, 0.55, 0.17, 0.17, 0.17, 0.59, 0.13, 0.15, 0.18, 0.13, 0.17, 0.13, 0.19, 0.12, 0.17, 0.19, 0.31, 0.29, 0.29, 0.29, 0.18, 0.18, 0.05, 0.18, 0.19, 0.15, 0.22, 0.13, 0.17]
   transrisk_static_work_sd::Array{Float64,1} = [transrisk_static_work_mean[i]*(transrisk_household_sd_average/transrisk_household_mean_average) for i=1:length(transrisk_static_work_mean)]
   transrisk_dynamic_work_mean::Array{Float64,1} = (transrisk_household_mean_average/0.5).*[0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.06, 0.06, 0.06, 0.33, 0.33, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.14, 0.21, 0.20, 0.14, 0.21, 0.25, 0.04, 0.35, 0.25, 0.25, 0.35, 0.33, 0.39, 0.39, 0.39, 0.36, 0.36, 0.25, 0.36, 0.35, 0.49, 0.56, 0.50, 0.25]
   transrisk_dynamic_work_sd::Array{Float64,1} = [transrisk_dynamic_work_mean[i]*(transrisk_household_sd_average/transrisk_household_mean_average) for i=1:length(transrisk_dynamic_work_mean)]

   # Baseline risk for social transmission, set relative to household transmission
   transrisk_social_mean::Float64 = (transrisk_household_mean_average/0.5)*0.355
   transrisk_social_sd::Float64 = transrisk_social_mean*(transrisk_household_sd_average/transrisk_household_mean_average)

   # Baseline risk for random transmission, set relative to household transmission
   transrisk_random_mean::Float64 = (transrisk_household_mean_average/0.5)*0.355
   transrisk_random_sd::Float64 = transrisk_random_mean*(transrisk_household_sd_average/transrisk_household_mean_average)

   # For workplaces that are covid-secure, scale the transmission risk
   # Zero would mean complete removal of transmission risk
   # One would mean no effect on the transmission risk
   CS_scale_transrisk::Array{Float64,1} = 0.2*ones(length(transrisk_static_work_mean))

   # probability of being asymptomatic
   probasymp_dist::Uniform{Float64} = Uniform(0.5,0.8)

   # Flag variable indicating if self-isolation is an active measure
   isolation::Int64 = 1
   symp_isoltime::Int64 = 10
   household_isoltime::Int64 = 14

   # Set proportion of population who will adhere
   adherence::Float64 = 0.7

   # Average social contacts on a workday & nonworkday
   n_social_mean_workday::Int64 = 1
   n_social_mean_nonworkday::Int64 = 5

   # Distribution of incubation period
   d_incub::Distribution = Erlang(6,0.88)

   # Distribution of infectiousness
   dist_infectivity::Array{Float64,1} =  [0.0369, 0.0491, 0.0835, 0.1190, 0.1439, 0.1497, 0.1354, 0.1076, 0.0757, 0.0476, 0.0269, 0.0138, 0.0064, 0.0044] # Corrected He (4 day pre-symptomatic period) & 10 day symptomatic period
   # dist_infectivity::Array{Float64,1} = [0.0379, 0.0504, 0.0857, 0.1220, 0.1475, 0.1535, 0.1388, 0.1103, 0.0776, 0.0488, 0.0276] # Corrected He (4 day pre-symptomatic period)
   # dist_infectivity::Array{Float64,1} = [0.1847,0.2371,0.1968,0.1396,0.0911,0.0566,0.0339,0.0199,0.0114]/sum([0.1847,0.2371,0.1968,0.1396,0.0911,0.0566,0.0339,0.0199,0.0114])

   # Distribution of delay in reporting symptomns (for those that do/eventually adhere)
   # Short delays. E.g. Symptom onset while at work & isolate next day = 1 day delay
   delay_adherence_pmf::Array{Float64,1} = [1.,0.,0.] # Note first entry is 0 day delay,
                                                         # second entry 1 day delay etc

   # Distribution of delay in household infection (for those that will be infected)
   # First entry corresponds to 0 days, second entry to 1 day, etc.
   #delay_household_infection_pmf::Array{Float64,1} = [1.,0.,0.,0.,0.,0.]
   delay_household_infection_pmf::Array{Float64,1} = dist_infectivity./sum(dist_infectivity)

   # At beginning of simulation, the proportion of student population that has
   # been infected previously
   recov_propn = 0.1
end

@with_kw mutable struct network_params
# used for parameters relating to the network generation

   # Number of nodes in the system
   n_nodes::Int64 = 0

   # Worker information array. Entry per worker.
   worker_nodes::Array{worker_params,1} = Array{worker_params,1}(undef,0)

   # Vector of vectors for workplace sizes
   # Vector per sector/worktype
   workplace_sizes::Array{Array{Int64,1},1} = Array{Int64,1}[]

   # Method used to generate contacts (can be "configuration" or "ER")
   network_generation_method::String = "configuration"

   # Method used to generate dynamic social contacts (can be "configuration", "ER" or "cluster")
   network_generation_method_dynamic_social::String = "cluster"

   # Connection probabilities with any other worker & socially
   # Consistent for everyone
   prob_anyworker_contact::Float64 = 0.0001*(10000/cmax)
   prob_social_contact::Float64 = 0.002*(10000/cmax)
   prob_random_contact::Float64 = 0.0001*(10000/cmax)

   # Household size distribution
   household_size_distribution::Array{Float64,1} = [0.314, 0.423, 0.075, 0.025, 0.006, 0.002]/sum([0.314, 0.423, 0.075, 0.025, 0.006, 0.002])

   # Sector specific worker contact probability
   prob_workertype_contact::Array{Float64,1} = Array{Float64,1}(undef,0)

   # Mean degree for workplaces in each sector
   dd_within_workplace::Array{Float64,1} = Array{Float64,1}(undef,0)

   # Degree distribution of contacts in different workplaces
   workplace_degree_distribution::Array{Distribution,1} = [Distributions.LogNormal(1.896,1.233)]

   # Probability of making contact with other workplace compared to within workplace
   between_workplace_contact_probs::Array{Float64,1} = [0.05]

   # Maximum number of dynamic contacts for an individual per day
   max_contacts_work_dynamic::Int64 = 100

   # Distribution of social group sizes
   social_group_size_distribution::Distribution = Distributions.LogNormal(3.14,1.41)
   # social_group_size_distribution::Distribution = Distributions.LogNormal(log(6),1.056)

   # Distribution of social contacts per day
   social_workday_dd::Distribution = Distributions.LogNormal(1.397,1.27)
   social_nonworkday_dd::Distribution = Distributions.LogNormal(1.536,1.153)

   # Maximum contacts allowed for social group
   max_contacts_social::Int64 = 100

   # Maximum group size for single meeting
   group_limit::Int64 = 100

   # Number of days to repeat same social contacts
   dynamic_time_frame::Int64 = 1

   n_groups_per_day_distribution::Distribution = Distributions.Poisson(1)

   # Probability of making contacts with f-o-f opposed to others
   friend_of_friend_prob::Float64 = 0.5

   # Whether or not to cluster daily social contacts
   cluster_social_contacts::Bool = true

   # mean number & standard deviation of dynamic contacts for each worker type
   dynamic_conts_mean::Array{Float64,1} = Array{Float64,1}(undef,0)
   dynamic_conts_sd::Array{Float64,1} = Array{Float64,1}(undef,0)

   # Degree distribution of dynamic contacts in different workplaces
   workplace_dynamic_degree_distribution::Array{Distribution,1} = [Distributions.LogNormal(1.262,1.315)]

   # Info on each workplace
   workplace_info::Array{Array{workplace_params,1},1} = Array{workplace_params,1}[]

   # Flag to indicate if sector is open (true) or closed (false)
   sector_open::Array{Bool,1} = Array{Bool,1}(undef,0)

   # Declares wether workplaces may have covid-secure (CS) status
   # If false, all workplaces are non-CS
   CS_active_flag::Bool = false

   # Specifies if workplace is covid-secure, max number of workplace contacts
   CS_team_size::Int64 = 2

   # Scaling on number of contacts for lockdowns
   social_contact_scaling::Array{Float64,1} = ones(Float64, cmax)
   random_contact_scaling::Array{Float64,1} = ones(Float64, cmax)

end

@with_kw mutable struct workplace_generation_params
# used for parameters relating to workplace generation
   workertypes::Int64 = 4
   workpercent::Array{Float64,1} = [1.,1.,1.,1.]
   workforce_proportion::Array{Float64,1} = [0.030, 0.127, 0.311, 0.532]
   workplace_size_mean::Array{Float64,1} = [6., 9., 11., 10.]
   workplace_size_sd::Array{Float64,1} = [14., 57., 122., 109. ]
   workplace_size_gen_fn::Function = workplace_size_sampled_from_empirical_pmf
end

@with_kw mutable struct sim_outputs
# used for outputs to be saved from the simulations
   endtime::Int64 = 0
   countfinal::Int64 = 0
   cmax::Int64 = 0
   # n_initial_asymp::Int64 = 0
   # n_initial_symp::Int64 = 0
   n_intervention_sets::Int64 = 0

   # 2D outputs
   numlat::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number latently infected
   numinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number infectious
   numrep::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number reporting infection
   prevlat::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for latently infected
   prevsymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for symptomatic infectious
   prevasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for asymptomatic infectious
   prevpresymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for recovereds
   newinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of new infections
   newasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of newly infected asymptomatics
   workersinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of workers that are newly infected
   workersasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of workers that are newly infected asymptomatics
   num_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating at each timepoint
   num_household_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating due to household members having symptoms at each timepoint
   num_symp_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating due to symptoms at each timepoint
   num_isolating_CTcause::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating as a contact of a positive test
   num_CT::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # total number of recallable contacts
   num_infected::Array{Int64,3} = zeros(Int64,cmax,countfinal,n_intervention_sets) # number of infections caused by that node
   dynamic_infection_count::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number infections over dynamic links
   Rt::Array{Float64,3} = zeros(Float64,endtime+1,countfinal,n_intervention_sets) # real time R value (number of secondary infections caused by nodes that were newly infected on that day)
   transmission_setting::Array{Int64,4} = zeros(Int64, endtime+1, countfinal, n_intervention_sets, 5) # records proportion of infections ocurring in each setting: social, work-static, work-dynamic, household, other
   num_init_infected::Array{Array{Int64,2}} = Array{Array{Int64,2}}(undef,countfinal) # number of infections caused by the initially infected nodes

   # 2D outputs. Testing related
   tests_performed::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets)

   # 3D outputs. Testing related
    # Slices: 1 - True positive; 2 - false negative; 3 - True negative; 4 - false positive.
   test_outcomes::Array{Int64,4} = zeros(Int64,endtime+1,countfinal,n_intervention_sets,4)

   # 1D outputs
   infected_by::Array{Int64,2} = zeros(Int64,cmax,n_intervention_sets)  # who each node was infected by
   var_num_infected::Array{Int64,2} = zeros(Int64,countfinal,n_intervention_sets) # variance in the number of infections caused by each node
   mean_init_generation_time::Array{Float64,2} = zeros(Float64,countfinal,n_intervention_sets) # mean initial generation time
end

@with_kw mutable struct intervention_data_feeds
   rep_inf_this_timestep::Array{Int64,1} = Array{Float64,1}[] #Indicator of whether a node reports symptomatic infection during current timestep
   numlat::Int64 = 0
   numinf::Int64 = 0
   numrep::Int64 = 0
   newinf::Int64 = 0
end

@with_kw mutable struct node_states
# used for the status of each node
   cmax::Int64 = 0
   timelat::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in latent state
   timeinf::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in infectious state
   timesymp::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in symptomatic state
   lattime::Array{Int64,1} = zeros(Int64,cmax) # time to spend in latent state
   inftime::Int64 = 4 # time to spend in infectious, pre-symptomatic state
   symptime::Int64 = 10 # time to spend in symptomatic state
   asymp::Array{Int64,1} = zeros(Int64,cmax) # whether the node will be asymptomatic
   rep_inf_this_timestep::Array{Int64,1} = zeros(Int64,cmax) # whether the node reports symptoms this timestep
   atwork::Array{Int64,2} = zeros(Int64,cmax,endtime) # whether the node is at work on each day
   hh_isolation::Array{Int64,1} = zeros(Int64,cmax) # Whether individual adheres to isolation guidance. (1) Yes. (0) No.
   delay_adherence::Array{Int64,1} = zeros(Int64,cmax) # Individial may not report symptoms immediately.
   acquired_infection::Array{Int64,1} = zeros(Int64,cmax) # time node acquired infection

   # Arrays to track isolation history
   hh_in_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Household isolation. Daily record for each student
   symp_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Symptomatic isolation. Daily record for each student
   CT_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Isolation by contact tracing. Daily record for each student
end

@with_kw mutable struct intervention_struct
   effects::Array{String,1} = ["none"]
   start_time::Int64 = 0
   sameday::Int64 = 3
   ton::Int64 = 1
   toff::Int64 = 0
   workpercent::Array{Float64} = Array{Float64}(undef)
   CT_parameters::CT_params = CT_params()
   scaling_work_static::Float64 = 1.
   scaling_work_dynamic::Float64 = 1.
   scaling_social::Float64 = 1.
   scaling_random::Float64 = 1.
   adherence::Float64 = 0.7
   social_contact_scaling::Float64 = 1.
   random_contact_scaling::Float64 = 1.
   CS_scale_transrisk_scalar::Float64 = 1.
   CS_team_size_intervention::Float64 = 1.
end
