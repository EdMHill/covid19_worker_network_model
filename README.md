# covid19_worker_network_model

This repository contains files for performing computational simulations of a network model framework to explore SARS-CoV-2 transmission amongst the workforce.

The code was developed for the analysis presented in the scientific paper "A network modelling approach to assess nonpharmaceutical disease controls in a worker population: An application to SARS-CoV-2" by Edward M. Hill, Benjamin D. Atkins, Matt J. Keeling, Louise Dyson and Michael J. Tildesley.

Publication details: Hill et al. (2021) A network modelling approach to assess non-pharmaceutical disease controls in a worker population: An application to SARS-CoV-2. *PLOS Computational Biology* 17(6): e1009058. doi: 10.1371/journal.pcbi.1009058. URL: https://doi.org/10.1371/journal.pcbi.1009058.

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://julialang.github.io/Pkg.jl/v1.5/environments/#Using-someone-else's-project-1

Please find below an explainer of the directory structure within this repository.

## Data
Directory containing workplace size data

## Results
Directory to store simulation outputs and plot scripts.

## src

**worker_pattern_model.jl**  
The main model run file. Pass requested configuration to worker_pattern_network_run (in include_files_network_model/main_function.jl) and saves outputs into an MAT file.

**include_files_network_model**  
Houses function files to be used when running the uni model.

- **main_function.jl**
    Outline of the code structure:  
    * Unpack required variables
    * Set the random number generator
    * Initialise workers into sectors, split by workplaces (supporting functions in network_generation_fns.jl)
    * Generate contacts within workplaces & households (supporting functions in network_generation_fns.jl)
    * Generate social contacts (workdays and non-workdays) (supporting functions in network_generation_fns.jl)
    * Generate workplace dynamic contacts (supporting functions in network_generation_fns.jl)
    * Generate random (social) dynamic contacts (supporting functions in network_generation_fns.jl)
    * Generate transmission risks (supporting functions in additional_fns.jl)
    * Initialise storage arrays and variables
    * Iterate over different intervention sets and run replicates:
        - Reinitialisation phase (supporting functions in additional_fns.jl)
        - Set course of infection times (supporting functions in additional_fns.jl)
        - Seed non-susceptible disease states (supporting functions in seed_initial_states_fn.jl)
        - Update output time series with initial conditions
        - Reset contact tracing variables (supporting functions in additional_fns.jl)
        - Iterate over time
            - Reinitialise variables at start of timestep
            - Implement any interventions
            - Assign outputs
            - Increment counters (supporting functions in additional_fns.jl)
            - Increment infection process (if come to the end of latent time, move to infectious time etc.). Includes household infection loop. (supporting functions in additional_fns.jl)
            - Record whether individuals are in isolation. Set workplace attendance status for timestep.
            - Transmit infections (supporting functions in additional_fns.jl)
            - Perform contact tracing (supporting functions in contact_tracing_fns.jl)
            - Assign prevalence & isolation outputs
            - Reactive workplace closure check
            - Run interventions (supporting functions in intervention_condition_affect_fns.jl)
    * Return outputs

- **additional_fns.jl**   
    Stash of functions to run the worker network model. Includes:  
    * populate arrays specifying when individuals are at the workplace (populate_atwork!)
    * time in state increment fn (increment_counters!)
    * load configs for sensitivity runs, generates (load_configs)
    * find_network_parameters (Load relevant network params based on number of work sectors requested)
    * transmission functions
    * functions to set up transmission rates within household for each individual
    * functions to reinitialise states at start of each simulation replicate
    * miscellaneous functions

- **seed_initial_states_fn.jl**   
    Function specifying how the initial disease states will be assigned each simulation replicate.

- **contact_tracing_fns.jl**  
    Functions that are used with the worker network model for performing contact tracing. Includes:
    * Lookup if homeday/workday contacts were active in contact tracing active window (lookup_homeday_workday_CT)
    * Get portion of dynamic contacts to be recallable (recallable_dynamic_contacts)
    * Check contacts made at work (get_worker_contacts)
    * Perform forward CT from an identified infector (forwardCT_from_infector! and trace_node!)

- **intervention_condition_affect_fns.jl**  
    Functions to implement/rescind trigger based interventions

- **network_generation_fns.jl**  
    Produces the network layers

- **parametertypes.jl**  
    Defines the parameter types used with the network model for workers. Fields accessible with dot notation. Example using type worker_params, with a variable named worker_nodes. worker_nodes.sector_ID accesses the value in the sector_ID field.    

- **workplace_size_generation_fns.jl**  
    Functions to generate workplace sizes for a given number of work sectors
