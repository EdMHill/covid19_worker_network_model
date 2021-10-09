#=
Purpose:
Stash functions that are used with the network model for modelling contact tracing
#-------------------------------------------------------------------------------
=#

#-------------------------------------------------------------------------------
# LOOKUP IF HOMEDAY/WORKDAY CONTACTS WERE ACTIVE IN CONTACT TRACING ACTIVE WINDOW
#-------------------------------------------------------------------------------

"""
    recallable_dynamic_contacts(worker_ID::Int64,
                                relevant_prev_days_for_CT::Int64,
                                atwork::Array{Int64,2},
                                time::Int64)

Find relevant workday & homeday contacts in time 'd' days before symptoms up to current day.

Inputs:
- `worker_ID`: Identification of worker
- relevant_prev_days_for_CT: Amount of days to be looked back over for CT purposes
- atwork: Binary array. Row per worker with entries (0) homeday, (1) workday
- time: Current timestep of the simulation

Outputs:
- homeday_flag::Bool - True if worker had a homeday at any point during relevant_prev_days_for_CT
- workday_flag::Bool- True if worker had a workday at any point during relevant_prev_days_for_CT

Location: contact\\_tracing\\_fns.jl
"""
function lookup_homeday_workday_CT(worker_ID::Int64,
                                    relevant_prev_days_for_CT::Int64,
                                    atwork::Array{Int64,2},
                                    time::Int64)
    # Take in contacts and worker schedule. Check if relevant


    prev_day_increment = 1 # Initialise counter to access previous days in atwork array
    homeday_flag = false; workday_flag = false;
    check_flag = true # Initialise check_flag, states whether you stay in loop
    while (check_flag == true) && (prev_day_increment <= relevant_prev_days_for_CT) && ((time - prev_day_increment) > 0)


        # Look back over relevant_prev_days_for_CT
        # Get if a workday (1), or home day (0)
        day_type_val = atwork[worker_ID,time-prev_day_increment]

        if day_type_val == 0
            homeday_flag = true
        else
            workday_flag = true
        end

        # Update previous day increment
        prev_day_increment += 1

        # Both workday and homeday satisfied
        if (homeday_flag == true) && (workday_flag == true)
            check_flag = false
        end
    end

    return homeday_flag::Bool,
            workday_flag::Bool
end

"""
    recallable_dynamic_contacts(args)

For a single specified node, collect expected dynamic contacts for a specified day and check if they occurred and are recallable.

Inputs: `worker_ID` - ID of node being contact traced,
        `dynamic_contact_record` - IDs of dynamic contacts for all nodes on all days,
        `dynamic_contacts_recalled_propn` - array defining proportion of contacts recalled x days ago,
        `daily_record_inisol` - indicates isolation status of each node each day,
        `daily_record_atworkplace` - indicates at workplace status of each node each day
        `time_to_check` - time that contacts to be checked took place,
        `prev_day_val` - number of days prior to current time that time_to_check is,
        `rng` - random number generator \n
Outputs: `recallable_contacts` - array of node IDs that have been recalled as contacts \n
Location: contact\\_tracing\\_fns.jl
"""
function recallable_dynamic_contacts(worker_ID::Int64,
                                        dynamic_contact_record::Array{Array{Int64,1},2},
                                        dynamic_contacts_recalled_propn::Array{Float64,1},
                                        daily_record_inisol::Array{Int64,2},
                                        daily_record_atworkplace::Array{Int64,2},
                                        time_to_check::Int64,
                                        prev_day_val::Int64,
                                        rng::MersenneTwister
                                        )

    # Get proportion of dynamic contacts from prev_day_val days ago to be retained
    if prev_day_val > length(dynamic_contacts_recalled_propn)
        # No contacts can possible be retained.
        # Return an empty vector
        recallable_contacts = Int64[]
    else
        # Get proportion of contacts expected to be retained from input distribution
        threshold = dynamic_contacts_recalled_propn[prev_day_val]

        # Get IDs of those dynamic contacts on required day
        all_dynamic_contacts = dynamic_contact_record[time_to_check,worker_ID]

        # Check if any dynamic contacts were isolating.
        # If so, will not be a recallable contact
        n_possible_dynamic_contacts = length(all_dynamic_contacts)
        recallable_contact_check = zeros(Int64,n_possible_dynamic_contacts)
        for contact_itr = 1:n_possible_dynamic_contacts
            # If not isolating, then contact did occur
            contact_ID = all_dynamic_contacts[contact_itr]
            if daily_record_inisol[time_to_check,contact_ID]==false
                # Contact occurred
                # Check if contact will be remembered
                r = rand(rng)
                if r<threshold
                    recallable_contact_check[contact_itr] = 1
                end
            end
        end

        # Construct vector of IDs of contacts that did actually occur
        n_recallable_contacts = sum(recallable_contact_check)
        recallable_contacts = zeros(Int64,n_recallable_contacts) # Initialise vector to store IDs of recallable contacts
        recallable_contact_itr = 1 # Initialise counter for vector assignment index
        for contact_itr = 1:n_possible_dynamic_contacts
            if recallable_contact_check[contact_itr] == true
                contact_ID = all_dynamic_contacts[contact_itr]
                recallable_contacts[recallable_contact_itr] = contact_ID
                recallable_contact_itr += 1 # Increment vector assignment index
            end
        end
    end

    return recallable_contacts::Array{Int64,1}
end

"""
    get_worker_contacts(args)

For a single specified node, go over usual workday contacts and check if they actually occurred on a given day.

Inputs: `worker_ID` - ID of node being contact traced,
        `possible_worker_contacts` - IDs of usual workplace contacts, to be checked,
        `daily_record_inisol` - indicates isolation status of each node each day,
        `daily_record_atworkplace` - indicates at workplace status of each node each day
        `time_to_check` - time that contacts to be checked took place,
        `prev_day_val` - number of days prior to current time that time_to_check is,
        `rng` - random number generator,
        `network_params` - NetworkParameters structure,
        `other_workplace_flag` - boolean flagging if tracing contacts in same or different workplace \n
Outputs: `workplace_contacts` - array of node IDs that have been confirmed as contacts \n
Location: contact\\_tracing\\_fns.jl
"""
function get_worker_contacts(worker_ID::Int64,
                                possible_worker_contacts::Array{Int64,1},
                                daily_record_inisol::Array{Int64,2},
                                daily_record_atworkplace::Array{Int64,2},
                                time_to_check::Int64,
                                prev_day_val::Int64,
                                rng::MersenneTwister,
                                network_parameters::network_params;
                                other_workplace_flag::Bool = false
                                )

    # Unpack parameters required
    if other_workplace_flag == true
        @unpack worker_nodes, workplace_info = network_parameters
    end

    # Check if any contacts were isolating and/or not at workplace that day
    # If so, will not be a contact
    n_possible_worker_contacts = length(possible_worker_contacts)
    contact_occur_check = zeros(Int64,n_possible_worker_contacts)
    for contact_itr = 1:n_possible_worker_contacts

        # If contact not isolating and at workplace, then contact can occur in same workplace.
        # If contact in another workplace, also need to check CS status
        contact_ID = possible_worker_contacts[contact_itr]
        if (daily_record_inisol[time_to_check,contact_ID]==false) &&
            (daily_record_atworkplace[time_to_check,contact_ID]==true)

            # Check if dealing with same workplace contact
            # or contact with other workplaces
            if other_workplace_flag == true
                # Contact in other workplace. Need to check CS status of that workplace

                # Get information on the contacts workplace
                contact_sector_ID = worker_nodes[contact_ID].sector_ID
                contact_workplace_ID = worker_nodes[contact_ID].workplace_ID
                contact_workplace_info = workplace_info[contact_sector_ID][contact_workplace_ID]

                # Get contacts workplace CS status
                # If not a CS setting, contact has taken place
                contact_workplace_CS_bool = contact_workplace_info.covid_secure
                if (contact_workplace_CS_bool==false)
                    contact_occur_check[contact_itr] = 1
                end
            else
                # Contact in same workplace
                contact_occur_check[contact_itr] = 1
            end
        end
    end

    # Construct vector of IDs of contacts that did actually occur
    n_occur_worker_contacts = sum(contact_occur_check)
    workplace_contacts = zeros(Int64,n_occur_worker_contacts) # Initialise vector to store IDs of recallable contacts
    contact_occur_itr = 1 # Initialise counter for vector assignment index
    for contact_itr = 1:n_possible_worker_contacts
        if contact_occur_check[contact_itr] == true
            contact_ID = possible_worker_contacts[contact_itr]
            workplace_contacts[contact_occur_itr] = contact_ID
            contact_occur_itr += 1 # Increment vector assignment index
        end
    end

    return workplace_contacts::Array{Int64,1}
end

"""
    forwardCT_from_infector!(args)

Perform forward contact-tracing from an identified infector, if infector reports symptoms and engages with contact-tracing.

Inputs: `infector_ID` - ID of infector of contact traced node,
        `... parameter structures ...`,
        `engage_with_CT` - array indicating whether each node engages with contact-tracing,
        `atwork` - work schedule for each node,
        `time` - current time in simulation,
        `count` - number of current replicate,
        `rng` - random number generator \n
Outputs: None \n
Location: contact\\_tracing\\_fns.jl
"""
function forwardCT_from_infector!(infector_ID::Int64,
                                CT_vars::contact_tracing_vars,
                                contacts::contacts_struct,
                                CT_parameters::CT_params,
                                engage_with_CT::Array{Bool,1},
                                atwork::Array{Int64,2},
                                time::Int64,
                                count::Int64,
                                hh_isolation::Array{Int64,1},
                                network_parameters::network_params,
                                rng::MersenneTwister,
                                infector_trace_count::Array{Int64,1})

@unpack dynamic_contacts_recalled_propn, social_contacts_recalled_propn, infector_engage_with_CT_prob = CT_parameters
@unpack Inds_to_be_contacted, Test_result_false_negative = CT_vars

    # Already contact traced? If not, do it now
    if (isassigned(Inds_to_be_contacted,infector_ID)) # Checks if infector_ID reported infection.
                                                      # If not, no forward contact tracing will be done from infector
        if !(length(Inds_to_be_contacted[infector_ID])>0) # Inds_to_be_contacted[infector_ID] being empty signifies
                                                         # infector reported symptoms, but returned false neg and/or
                                                         # did not engage in contact tracing
            # check if the infector will engage with CT
            if (((engage_with_CT[infector_ID] == false) && (rand(rng)<infector_engage_with_CT_prob))
                || (engage_with_CT[infector_ID] == true)) &&
                (Test_result_false_negative[infector_ID] == false)
                        # didn't engage before but does now or already willing to engage
                        # and then checks infector has not previously tested negative

                # println("Infector engages")
                infector_trace_count[1] += 1

                trace_node!(infector_ID,time,CT_vars,contacts,CT_parameters,network_parameters,rng)
            end
        end
    end
end

"""
    trace_node!(node_itr::Int64,time::Int64,CT_vars::contact_tracing_vars,
                contacts::contacts_struct,CT_parameters::CT_params,
                network_parameters::network_params,rng)

Perform contact tracing from an index node.

Outputs: None \n
Location: contact\\_tracing\\_fns.jl
"""
function trace_node!(node_itr::Int64,time::Int64,CT_vars::contact_tracing_vars,
    contacts::contacts_struct,
    CT_parameters::CT_params,network_parameters::network_params,rng)
@unpack worker_nodes, workplace_info, CS_active_flag = network_parameters

    # Preload workplace contacts and other worker contacts
    if CS_active_flag == true
        possible_same_workplace_contacts = contacts.work_contacts_same_workplace_CS[node_itr]
    else
        possible_same_workplace_contacts = contacts.work_contacts_same_workplace[node_itr]
    end
    possible_other_workplace_contacts = contacts.work_contacts_other_workplace[node_itr]

    # Preload workplace info
    workertype_ID = worker_nodes[node_itr].sector_ID
    workplace_ID = worker_nodes[node_itr].workplace_ID
    current_workplace_info = workplace_info[workertype_ID][workplace_ID]

    # Get contacts that will be contacted
    for time_itr = 1:CT_vars.relevant_prev_days_for_CT[node_itr]
        if time-time_itr > 0 # Can't look back before the simulation started

            # Get previous time being checked
            time_to_check = time-time_itr

            # If at workplace, get contacts made that day
            if (contacts.daily_record_atworkplace[time_to_check,node_itr]==1)

                # Note, function "get_worker_contacts"
                # resides in contact_tracing_fns.jl
                # Workplace contacts
                workplace_contacts = get_worker_contacts(node_itr,
                                                            possible_same_workplace_contacts,
                                                            contacts.daily_record_inisol,
                                                            contacts.daily_record_atworkplace,
                                                            time_to_check,
                                                            time_itr,
                                                            rng,
                                                            network_parameters)

                # If there are workplace contacts,
                # add to vector tracking traceable contacts
                if !isempty(workplace_contacts)
                    append!(CT_vars.Inds_to_be_contacted[node_itr],workplace_contacts)
                end
                # Workers in a CS settings will not make any non-workplace
                # worker contacts
                if (current_workplace_info.covid_secure == false)
                    # Worker contacts at other workplaces
                    other_workplace_flag_val = true
                    other_workplace_contacts = get_worker_contacts(node_itr,
                                                                possible_other_workplace_contacts,
                                                                contacts.daily_record_inisol,
                                                                contacts.daily_record_atworkplace,
                                                                time_to_check,
                                                                time_itr,
                                                                rng,
                                                                network_parameters,
                                                                other_workplace_flag=other_workplace_flag_val
                                                                )

                    # If there are workday contacts with
                    # inds at other workplaces,
                    # add to vector tracking traceable contacts
                    if !isempty(other_workplace_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],other_workplace_contacts)
                    end
                end

                # Dynamic contacts as part of job role at workplace
                # Check at workplace on current timestep
                atwork_dynamic_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                    contacts.dynamic_worker_contacts,
                                                                    CT_parameters.dynamic_contacts_recalled_propn,
                                                                    contacts.daily_record_inisol,
                                                                    contacts.daily_record_atworkplace,
                                                                    time_to_check,
                                                                    time_itr,
                                                                    rng) # in contact_tracing_fns.jl

                # If there are recallable dynamic contacts, add to vector tracking traceable contacts
                if !isempty(atwork_dynamic_recallable_contacts)
                    append!(CT_vars.Inds_to_be_contacted[node_itr],atwork_dynamic_recallable_contacts)
                end
            end

            # Social contacts check
            # Check node not in isolation, and has
            # potentially made social contacts
            if (contacts.daily_record_inisol[time_to_check,node_itr]==false) &&
                (contacts.social_contacts_per_node[node_itr] > 0)

                # Workday social contacts
                if (contacts.daily_record_atworkplace[time_to_check,node_itr]==true)
                    workday_social_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                        contacts.workday_social_contacts_by_day,
                                                                        CT_parameters.social_contacts_recalled_propn,
                                                                        contacts.daily_record_inisol,
                                                                        contacts.daily_record_atworkplace,
                                                                        time_to_check,
                                                                        time_itr,
                                                                        rng) # in contact_tracing_fns.jl

                    if !isempty(workday_social_recallable_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],workday_social_recallable_contacts)
                    end

                else
                    # Non-workday social contacts
                    nonworkday_social_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                    contacts.nonworkday_social_contacts_by_day,
                                                                    CT_parameters.social_contacts_recalled_propn,
                                                                    contacts.daily_record_inisol,
                                                                    contacts.daily_record_atworkplace,
                                                                    time_to_check,
                                                                    time_itr,
                                                                    rng) # in contact_tracing_fns.jl

                    # If there are recallable dynamic contacts, add to vector tracking traceable contacts
                    if !isempty(nonworkday_social_recallable_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],nonworkday_social_recallable_contacts)
                    end
                end
            end
        end
    end
end
