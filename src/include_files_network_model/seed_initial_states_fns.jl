#=
Purpose:
Store functions used to seed nodes in non-susceptible states at beginning of
simulation replicates

Fns to select nodes to begin in given non-susceptible state
- choose_from_all_popn  (select from the entire population)

Fns to get number of nodes to be seeded in each non-susceptible state
- set_ten_initial_infected
- seed_states_with_uncertainty
#-------------------------------------------------------------------------------
=#

#-------------------------------------------------------------------------------
# FUNCTIONS TO SELECT NODES TO BEGIN IN GIVEN NON-SUSCEPTIBLE STATE
#-------------------------------------------------------------------------------

# select from the entire population
function choose_from_all_popn(rng::MersenneTwister,
                                    n_nodes::Int64,
                                    states::node_states,
                                    probasymp::Float64,
                                    infected_by::Array{Int64,1},
                                    n_initial_latent::Int64,
                                    n_initial_asymp::Int64,
                                    n_initial_symp::Int64,
                                    n_initial_rec::Int64,
                                    initialise_start_disease_state_flag::Bool)
# Inputs:
# rng::MersenneTwister - The random number generator
# n_nodes::Int64 - Numer of nodes in the network
# states::node_states - Record of info per node
# probasymp::Float64 - Probability of infected case being asymptomatic
# infected_by::Array{Int64,1} - Record of who each node was infected by
# n_initial_latent, n_initial_asymp, n_initial_symp, n_initial_rec
#   -  Numbers to be seeded in the named disease state
# initialise_start_disease_state_flag::Bool - If true, those assigned to non-susceptible state are initialised
#                   as having just entered that disease state (no elapsed time prior occurred)

# Ouputs:
# Directly mutates inputs states and infected_by


    # Initialise initial asymp infectious nodes
    for seed_inf_itr = 1:n_initial_asymp
       chosen_node = ceil(Int64,rand(rng)*n_nodes)

       valid_asymp_node = false
       while valid_asymp_node == false

           if (states.timeinf[chosen_node] > 0) ||
                (states.timesymp[chosen_node] > 0) # Chosen node already selected as an asymptomatic
              # Resample
              chosen_node = ceil(Int64,rand(rng)*n_nodes)
           else
              # Chosen node has asymptomatic status.
              valid_asymp_node = true
           end
       end

       # Set node as an asymptomatic
       states.asymp[chosen_node] = 1

       # Update states for chosen node
       states.timelat[chosen_node] = -1
       infected_by[chosen_node] = -1

       # Update time elapsed in state
       if initialise_start_disease_state_flag == true
           # Just entered infectious state
           states.timeinf[chosen_node] = 1
           states.acquired_infection[chosen_node] = -states.lattime[chosen_node]
       else
           # Sample time already elapsed as infected
           infected_period_length = states.inftime + states.symptime
           infected_time_elapsed::Int64 = rand(rng,1:infected_period_length)

           if infected_time_elapsed > states.inftime  # Elapsed infection time beyond pre-symptomatic stage
               states.timesymp[chosen_node] = infected_time_elapsed - states.inftime
               states.timeinf[chosen_node] = -1
           else
               states.timeinf[chosen_node] = infected_time_elapsed
           end
           states.acquired_infection[chosen_node] = -states.lattime[chosen_node] - (infected_time_elapsed - 1)
       end
    end

    # Initialise initial symptomatic infectious nodes
    for seed_inf_itr = 1:n_initial_symp
            chosen_node = ceil(Int64,rand(rng)*n_nodes)

            valid_symp_node = false
            while valid_symp_node == false
                if (states.asymp[chosen_node] == 1) ||  # Chosen node does not have symptomatic status.
                        (states.timeinf[chosen_node] > 0) ||
                        (states.timesymp[chosen_node] > 0) # Chosen node already selected as infectious
                   # Resample
                   chosen_node = ceil(Int64,rand(rng)*n_nodes)
                else
                   # Chosen node has symptomatic status.
                   valid_symp_node = true
                end
            end

            # Update states for chosen node
            states.timelat[chosen_node] = -1
            infected_by[chosen_node] = -1

            # Update time elapsed in state
            if initialise_start_disease_state_flag == true
                # Just entered infectious state
                states.timeinf[chosen_node] = 1
                states.acquired_infection[chosen_node] = -states.lattime[chosen_node]
            else
                # Sample time already elapsed as infected (from pre-symptomatic period)
                presymp_time_elapsed::Int64 = rand(rng,1:states.inftime)

                states.timeinf[chosen_node] = presymp_time_elapsed
                states.acquired_infection[chosen_node] = -states.lattime[chosen_node] - (presymp_time_elapsed - 1)
             end
    end

    # Initialise initial latent nodes
    for seed_latent_itr = 1:n_initial_latent
       chosen_latent_node = ceil(Int64,rand(rng)*n_nodes)

       # Check if node already set to be an initial infected or initial latent
       valid_latent_node = false
       while valid_latent_node == false
              if (states.timeinf[chosen_latent_node] > 0) ||
                  (states.timesymp[chosen_latent_node] > 0) || # Initial infected condition
                   (states.timelat[chosen_latent_node] > 0)  # Initial latent condition

                   # Redraw sample if already set to be an initial infected or initial recovered
                   chosen_latent_node = ceil(Int64,rand(rng)*n_nodes)
              else
                   # Valid node now drawn. Update flag
                   valid_latent_node = true
              end
       end

       # Update states variables
       infected_by[chosen_latent_node] = -2     # Unknown infector (someone outside popn). Set to -2 (rather than -1) to differentiate from initial infectious nodes.

       # Update time elapsed in state
       if initialise_start_disease_state_flag == true
           # Just entered infectious state
           states.timelat[chosen_node] = 1
           states.acquired_infection[chosen_node] = 0
       else
           # Sample time spent latent from length of lattime
           latent_period_length = states.lattime[chosen_latent_node]
           latent_time_elapsed::Int64 = rand(rng,1:latent_period_length)

           # Update states variables
           states.timelat[chosen_latent_node] = latent_time_elapsed
           states.acquired_infection[chosen_latent_node] = -(latent_time_elapsed - 1)
       end

       # Check if infection will be asymptomatic
       if rand(rng) < probasymp
           states.asymp[chosen_latent_node] = 1
       end
    end

    # Initialise nodes that are already recovered
    for seed_rec_itr = 1:n_initial_rec
        chosen_rec_node = ceil(Int64,rand(rng)*n_nodes)

        # Check if node already set to be an initial infected, initial latent
        # or initial recovered
        valid_node = false
        while valid_node == false
            if (states.timeinf[chosen_rec_node] > 0) ||
                (states.timesymp[chosen_rec_node] > 0) || # Initial infected condition
                  (states.timelat[chosen_rec_node] > 0) || # Initial latent condition
                ( (states.timelat[chosen_rec_node] == -1) && # Initial recovered conditions
                    (states.timeinf[chosen_rec_node] == -1) &&
                    (states.timesymp[chosen_rec_node] == -1) )

                # Redraw sample if already set to be an initial infected or initial recovered
                chosen_rec_node = ceil(Int64,rand(rng)*n_nodes)
            else
                # Valid node now drawn. Update flag
                valid_node = true
            end
        end

        # Set infection state variables to -1 (as previously progressed through them)
        states.timelat[chosen_rec_node] = -1
        states.timeinf[chosen_rec_node] = -1
        states.timesymp[chosen_rec_node] = -1
        infected_by[chosen_rec_node] = -3      # Unknown infector (someone outside student popn).
    end

    return nothing
end

#-------------------------------------------------------------------------------
# FUNCTIONS TO GET NUMBER OF NODES SEEDED IN EACH NON-SUSCEPTIBLE STATE
#-------------------------------------------------------------------------------

# Standard function setup
# Inputs
# rng::MersenneTwister - The random number generator
# n_nodes::Int64 - Numer of nodes in the network
# states::node_states - Record of info per node
# infected_by::Array{Int64,1} - Record of who each node was infected by
# recov_propn::Float64 - Proportion of population to begin in recovered state

# Outputs
# n_initial_latent, n_initial_asymp, n_initial_symp, n_initial_rec
#   -  Numbers to be seeded in the named disease state

# Set up with 9 asymps & 1 symp
function set_ten_initial_infected(rng::MersenneTwister,
                                    n_nodes::Int64,
                                    states::node_states,
                                    probasymp::Float64,
                                    infected_by::Array{Int64,1},
                                    recov_propn::Float64)

    # Set initial infected counts
    n_initial_latent = 0
    n_initial_asymp = 9
    n_initial_symp = 1

    # Set initial recovereds based on recov_propn
    n_initial_rec = ceil(Int64,recov_propn*n_nodes)

    # Select nodes from the population
    initialise_start_disease_state_flag = true
    choose_from_all_popn(rng,
                            n_nodes,
                            states,
                            probasymp,
                            infected_by,
                            n_initial_latent,
                            n_initial_asymp,
                            n_initial_symp,
                            n_initial_rec,
                            initialise_start_disease_state_flag)

    return n_initial_latent::Int64,
        n_initial_asymp::Int64,
        n_initial_symp::Int64,
        n_initial_rec::Int64
end

# Example to generate counts from distributions
function seed_states_with_uncertainty(rng::MersenneTwister,
                                    n_nodes::Int64,
                                    states::node_states,
                                    probasymp::Float64,
                                    infected_by::Array{Int64,1},
                                    recov_propn::Float64)

    # Set distributions to draw counts from
    d_asymp = Uniform(5,8)

    # Set initial infected counts
    n_initial_latent = 0
    n_initial_asymp = round(Int64,rand(rng,d_asymp))
    n_initial_symp = 10 - n_initial_asymp

    # Set initial recovereds based on recov_propn
    n_initial_rec = ceil(Int64,recov_propn*n_nodes)

    # Select nodes from the population
    initialise_start_disease_state_flag = true
    choose_from_all_popn(rng,
                            n_nodes,
                            states,
                            probasymp,
                            infected_by,
                            n_initial_latent,
                            n_initial_asymp,
                            n_initial_symp,
                            n_initial_rec,
                            initialise_start_disease_state_flag)

    return n_initial_latent,
        n_initial_asymp,
        n_initial_symp,
        n_initial_rec
end
