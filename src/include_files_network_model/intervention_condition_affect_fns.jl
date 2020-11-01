"""
Purpose:
Define conditions and affects of closure/reopening interventions
Motivated by a wish to assess impact of sector closures
"""


"""
Condition functions. Use time and/or disease state variable values
For example, if incidence of symptomatic infection exceeds given level, apply restriction
"""

function condition_close_example(intervention_trigger_input_data::intervention_data_feeds,
                                        time::Int64,
                                        network_parameters::network_params)
# Output:
#   output_bool - Denotes whether condition outcome was true or false.

   @unpack n_nodes = network_parameters
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition surpassed
   daily_incidence = sum(rep_inf_this_timestep)/n_nodes
   if (daily_incidence > 0.01)
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end

function condition_open_example(intervention_trigger_input_data::intervention_data_feeds,
                                        time::Int64,
                                        network_parameters::network_params)
 # Output:
#   output_bool - Denotes whether condition outcome was true or false.

   @unpack n_nodes  = network_parameters
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition declined to low enough level
   daily_incidence = sum(rep_inf_this_timestep)/n_nodes
   if daily_incidence < 0.001
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end



"""
Affect functions
For example, closing a particular sector
"""

# Close sectors
function affect_close_example!(network_parameters::network_params)

   @unpack workplace_info, sector_open = network_parameters

   # Get number of sectors in use
   n_sectors = length(workplace_info)

   # Throw an error if number of sectors is not compatible with this function
   if n_sectors != 41
      error("Number of sectors in use is $(n_sectors). Selected affect function requires use of 41 sectors.")
   end

   # The sectors to be closed
   close_sector_IDs = [9,39]
   n_sectors_closed = length(close_sector_IDs)

   # Iterate over each sector to be closed
   for close_sector_itr = 1:n_sectors_closed
      # Get sector of interest
      close_sector_ID = close_sector_IDs[close_sector_itr]

      if (sector_open[close_sector_ID] == true)
         # Get number of workplaces in that sector
         n_workplaces = length(workplace_info[close_sector_ID])

         # Set status of each workplace in that sector to be closed
         for workplace_itr = 1:n_workplaces
            workplace_info[close_sector_ID][workplace_itr].workplace_open = false
         end

         # Update sector_open field
         sector_open[close_sector_ID] = false
      end
   end
end

# Open sectors
function affect_open_example!(network_parameters::network_params)

   @unpack workplace_info, sector_open = network_parameters

   # Get number of sectors in use
   n_sectors = length(workplace_info)

   # Throw an error if number of sectors is not compatible with this function
   if n_sectors != 41
      error("Number of sectors in use is $(n_sectors). Selected affect function requires use of 41 sectors.")
   end

   # The sectors to be opened
   open_sector_IDs = [9,39]
   n_sectors_open = length(open_sector_IDs)

   # Iterate over each sector to be closed
   for open_sector_itr = 1:n_sectors_open
      # Get sector of interest and number of workplaces in that sector
      open_sector_ID = open_sector_IDs[open_sector_itr]

      if (sector_open[open_sector_ID] == false)
         # Get number of workplaces in that sector
         n_workplaces = length(workplace_info[open_sector_ID])

         # Set status of each workplace in that sector to be closed
         for workplace_itr = 1:n_workplaces
            workplace_info[open_sector_ID][workplace_itr].workplace_open = true
         end

         # Update sector_open field
         sector_open[open_sector_ID] = true
      end
   end
end


function affect_intervention!(intervention::intervention_struct,
                                 rng::MersenneTwister,
                                 RNGseed::Int64,
                                cmax::Int64,
                                ton::Int64,
                                toff::Int64,
                                infection_parameters::infection_params,
                                infection_parameters_preintervention::infection_params,
                                sameday::Int64,
                                countfinal::Int64,
                                endtime::Int64,
                                CT_parameters::CT_params,
                                CT_vars::contact_tracing_vars,
                                contacts::contacts_struct,
                                network_parameters::network_params,
                                nodes_by_workplace::Array{Array{Array{Int64,1},1},1},
                                workplace_generation_parameters::workplace_generation_params,
                                workplace_closure_active::Bool,
                                states::node_states,
                                workplace_thresholds::Array{Array{Int64,1},1},
                                workertypes::Int64,
                                workplace_sizes::Array{Array{Int64,1},1},
                                assign_workplace_static_transrisk_fn::Function,
                                assign_workplace_dynamic_transrisk_fn::Function,
                                assign_social_transrisk_fn::Function,
                                assign_random_transrisk_fn::Function)

      # Info required for intervention:
      #    - sameday
      #    - ton
      #    - toff
      #    - CT params
      #    - workpercent
      #    - timing
      #    - transscaling
      #    - CS params
      #    - infection params
      #    - network params
      #    - workplace generation params

      # Redefine who has returned to work
      if "workpercent" ∈ intervention.effects
         redefine_returned_to_work!(rng,
                                    cmax,
                                    network_parameters.worker_nodes,
                                    intervention,
                                    workplace_generation_parameters)
      end

      # Rescale transrisk
      if "transrisk" ∈ intervention.effects
         redefine_transmission_risks!(RNGseed,
                                       intervention,
                                       infection_parameters,
                                       infection_parameters_preintervention,
                                       network_parameters,
                                       workertypes,
                                       assign_workplace_static_transrisk_fn,
                                       assign_workplace_dynamic_transrisk_fn,
                                       assign_social_transrisk_fn,
                                       assign_random_transrisk_fn)
      end

      # Contact tracing variables
      if "contact_tracing" ∈ intervention.effects
         redefine_contact_tracing!(rng,
                                    cmax,
                                    intervention,
                                    CT_parameters,
                                    CT_vars,
                                    workplace_thresholds,
                                    workertypes,
                                    workplace_sizes)
      end

      # Adherence
      if "adherence" ∈ intervention.effects
         redefine_adherence!(rng,
                              cmax,
                              intervention,
                              states,
                              infection_parameters)

         for node_itr = 1:cmax
            if states.hh_isolation[node_itr] == 1
               network_parameters.social_contact_scaling[node_itr] = minimum(network_parameters.social_contact_scaling)
               network_parameters.random_contact_scaling[node_itr] = minimum(network_parameters.random_contact_scaling)
            else
               network_parameters.social_contact_scaling[node_itr] = 1.
               network_parameters.random_contact_scaling[node_itr] = 1.
            end
         end
      end

      # Contact structure
      if "contact_structure" ∈ intervention.effects
         for node_itr = 1:cmax
            if states.hh_isolation[node_itr] == 1
               network_parameters.social_contact_scaling[node_itr] = intervention.social_contact_scaling
               network_parameters.random_contact_scaling[node_itr] = intervention.random_contact_scaling
            else
               network_parameters.social_contact_scaling[node_itr] = 1.
               network_parameters.random_contact_scaling[node_itr] = 1.
            end
         end
      end

      # COVID-secure measures
      if "COVID_secure" ∈ intervention.effects
         redefine_COVID_secure!(rng,
                                 cmax,
                                 intervention,
                                 infection_parameters,
                                 workertypes,
                                 network_parameters,
                                 contacts,
                                 nodes_by_workplace
                                 )
      end

      return nothing
end


function redefine_COVID_secure!(rng::MersenneTwister,
                                 cmax::Int64,
                                 intervention::intervention_struct,
                                 infection_parameters::infection_params,
                                 workertypes::Int64,
                                 network_parameters::network_params,
                                 contacts::contacts_struct,
                                 nodes_by_workplace::Array{Array{Array{Int64,1},1},1}
                                 )
   # Redefine transmission risk across contacts in COVID-secure settings
   infection_parameters.CS_scale_transrisk = intervention.CS_scale_transrisk_scalar*ones(workertypes)

   # Reconstruct work contacts to conform to capped team sizes
   # Turn on COVID-secure measures
   network_parameters.CS_active_flag = true

   # Initialise CS related variables
   contacts.work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
   contacts.work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
   for node_itr = 1:cmax
      contacts.work_contacts_same_workplace_CS[node_itr] = Int64[]
   end

   # Construct COVID-secure work contacts
   network_parameters.CS_team_size = intervention.CS_team_size_intervention
   CS_workplace_generation!(network_parameters.worker_nodes,
                              nodes_by_workplace,
                              contacts.work_contacts_same_workplace_CS,
                              contacts.work_contacts_same_workplace_per_node_CS,
                              network_parameters,
                              rng)

   # Set all workplaces to be COVID-secure
   for worker_grp_idx = 1:workertypes
      workplace_count = length(network_parameters.workplace_info[worker_grp_idx])
      for workplace_itr = 1:workplace_count
         network_parameters.workplace_info[worker_grp_idx][workplace_itr].covid_secure = true
      end
   end

   return nothing
end

function redefine_returned_to_work!(rng::MersenneTwister,
                                    cmax::Int64,
                                    worker_nodes::Array{worker_params,1},
                                    intervention::intervention_struct,
                                    workplace_generation_parameters::workplace_generation_params)

   @unpack workertypes, workpercent = workplace_generation_parameters

   if length(intervention.workpercent) != workertypes
      error("Number of sectors specified in intervention is inconsistent.")
   end

   # Initialise array to store relative change in workpercent for each sector
   workpercent_change_prob::Array{Float64,1} = Array{Float64,1}(undef, workertypes)
   workpercent_change_direction::Array{Int64,1} = Array{Int64,1}(undef, workertypes)

   # Iterate over sectors
   for sector_itr = 1:workertypes
      workpercent_diff = intervention.workpercent[sector_itr] - workpercent[sector_itr]
      if workpercent_diff > 0
         workpercent_change_prob[sector_itr] = workpercent_diff / (1 - workpercent[sector_itr])
         workpercent_change_direction[sector_itr] = 1
      elseif workpercent_diff < 0
         workpercent_change_prob[sector_itr] = -1*workpercent_diff / workpercent[sector_itr]
         workpercent_change_direction[sector_itr] = -1
      else
         workpercent_change_prob[sector_itr] = 0
         workpercent_change_direction[sector_itr] = 0
      end
   end

   # Iterate over all nodes
   for node_itr = 1:cmax
      # find details of current node
      sector_id = worker_nodes[node_itr].sector_ID
      returned_to_work = worker_nodes[node_itr].returned_to_work

      # if increased workforce proportion, iterate over workers currently WFH and send some back to work
      if (workpercent_change_direction[sector_id] > 0) && (returned_to_work == 0)
         if rand(rng) < workpercent_change_prob[sector_id]
            worker_nodes[node_itr].returned_to_work = 1
         end
      # If decreased workforce proportion, do opposite
      elseif (workpercent_change_direction[sector_id] < 0) && (returned_to_work == 1)
         if rand(rng) < workpercent_change_prob[sector_id]
            worker_nodes[node_itr].returned_to_work = 0
         end
      end

   end

   # Update workplace_generation_parameters
   workplace_generation_parameters.workpercent = intervention.workpercent

   return nothing
end


function redefine_contact_tracing!(rng::MersenneTwister,
                                    cmax::Int64,
                                    intervention::intervention_struct,
                                    CT_parameters::CT_params,
                                    CT_vars::contact_tracing_vars,
                                    workplace_thresholds::Array{Array{Int64,1},1},
                                    workertypes::Int64,
                                    workplace_sizes::Array{Array{Int64,1},1})

   CT_parameters.contact_tracing_active = intervention.CT_parameters.contact_tracing_active
   CT_parameters.CT_engagement = intervention.CT_parameters.CT_engagement
   CT_parameters.prob_backwards_CT = intervention.CT_parameters.prob_backwards_CT
   CT_parameters.workplace_CT_threshold = intervention.CT_parameters.workplace_CT_threshold
   CT_parameters.perform_CT_from_infector = intervention.CT_parameters.perform_CT_from_infector
   CT_parameters.workplace_closure_active = intervention.CT_parameters.workplace_closure_active

   for ii = 1:cmax
      engage_with_CT_rand = rand(rng)
       if engage_with_CT_rand < CT_parameters.CT_engagement # engage with contact tracing
           CT_vars.Engage_with_CT[ii] = true
       else # do not engage with contact tracing
           CT_vars.Engage_with_CT[ii] = false
       end
    end

    for worktypeID = 1:workertypes
          for work_ID = 1:length(workplace_sizes[worktypeID])
             workplace_thresholds[worktypeID][work_ID] = ceil(Int64,workplace_sizes[worktypeID][work_ID]*CT_parameters.workplace_CT_threshold)
          end
     end


   return nothing
end


function redefine_transmission_risks!(RNGseed::Int64,
                                       intervention::intervention_struct,
                                       infection_parameters::infection_params,
                                       infection_parameters_preintervention::infection_params,
                                       network_parameters::network_params,
                                       workertypes::Int64,
                                       assign_workplace_static_transrisk_fn::Function,
                                       assign_workplace_dynamic_transrisk_fn::Function,
                                       assign_social_transrisk_fn::Function,
                                       assign_random_transrisk_fn::Function)

   infection_parameters.transrisk_static_work_mean = infection_parameters_preintervention.transrisk_static_work_mean*intervention.scaling_work_static
   infection_parameters.transrisk_dynamic_work_mean = infection_parameters_preintervention.transrisk_dynamic_work_mean*intervention.scaling_work_dynamic
   infection_parameters.transrisk_social_mean = infection_parameters_preintervention.transrisk_social_mean*intervention.scaling_social
   infection_parameters.transrisk_random_mean = infection_parameters_preintervention.transrisk_random_mean*intervention.scaling_random

   infection_parameters.transrisk_static_work_sd = infection_parameters_preintervention.transrisk_static_work_sd*intervention.scaling_work_static
   infection_parameters.transrisk_dynamic_work_sd = infection_parameters_preintervention.transrisk_dynamic_work_sd*intervention.scaling_work_dynamic
   infection_parameters.transrisk_social_sd = infection_parameters_preintervention.transrisk_social_sd*intervention.scaling_social
   infection_parameters.transrisk_random_sd = infection_parameters_preintervention.transrisk_random_sd*intervention.scaling_random

   assign_workplace_static_transrisk_fn(RNGseed,
                                        network_parameters,
                                        workertypes,
                                        infection_parameters.transrisk_static_work_mean,
                                        infection_parameters.transrisk_static_work_sd)

    assign_workplace_dynamic_transrisk_fn(RNGseed,
                                        network_parameters,
                                        workertypes,
                                        infection_parameters.transrisk_dynamic_work_mean,
                                        infection_parameters.transrisk_dynamic_work_sd)

    assign_social_transrisk_fn(RNGseed,
                                network_parameters,
                                infection_parameters.transrisk_social_mean,
                                infection_parameters.transrisk_social_sd)

    assign_random_transrisk_fn(RNGseed,
                                network_parameters,
                                infection_parameters.transrisk_random_mean,
                                infection_parameters.transrisk_random_sd)

   return nothing
end

function redefine_adherence!(rng::MersenneTwister,
                              cmax::Int64,
                              intervention::intervention_struct,
                              states::node_states,
                              infection_parameters::infection_params)

   current_adherence = infection_parameters.adherence
   new_adherence = intervention.adherence

   adherence_change_diff = new_adherence - current_adherence

   if adherence_change_diff > 0
      adherence_change_prob = adherence_change_diff / (1 - current_adherence)

      for node_itr = 1:cmax
         if (states.hh_isolation[node_itr] == 0) && (rand(rng) < adherence_change_prob)
            states.hh_isolation[node_itr] = 1
         end
      end

   elseif adherence_change_diff < 0
      adherence_change_prob = (-1*adherence_change_diff) / current_adherence

      for node_itr = 1:cmax
         if (states.hh_isolation[node_itr] == 1) && (rand(rng) < adherence_change_prob)
            states.hh_isolation[node_itr] = 0
         end
      end
   end

   return nothing
end
