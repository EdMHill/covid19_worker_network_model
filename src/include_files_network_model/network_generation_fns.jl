"""
Purpose:
Functions to produce the network layers

- generate_workplaces_and_allocate_workers (Give workplace to each worker)
- generate_contacts                        (Generate the networks)
- CS_workplace_generation!                 (Covid-secure workplace generation)
- configuration_model!                     (Network construction using configuration model)
- ER_model!                                (Network construction using erdos-renyi)
- generate_dynamic_worker_contacts         (Premake dynamic worker contacts,
                                            to be loaded in ahead of simulation)
- generate_social_contacts_each_day         (Premake social contacts (from friendship group) made each day,
                                            to be loaded in ahead of simulation)
- generate_random_contacts                  (premake social contacts outside friendship group)
"""

## Functions to produce the initial network layers ##

# Initialise vectors determining how many workers there are per node and which workplace each worker belongs to
function generate_workplaces_and_allocate_workers(cmax::Int64,
                                            workplace_generation_parameters::workplace_generation_params,
                                            RNGseed::Int64, rng::MersenneTwister, CS_active_flag::Bool)

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Get workplace params based on number of work sectors in use
    @unpack workertypes, workpercent, workforce_proportion, workplace_size_mean, workplace_size_sd = workplace_generation_parameters

    # total number of workers in each worker type
    worker_numbers = [round(Int64,cmax*workforce_proportion[i]) for i=1:length(workforce_proportion)]
    # correct for rounding
    if sum(worker_numbers) != cmax
        diff = sum(worker_numbers) - cmax
        worker_numbers[end] -= diff
    end

    # Initialise places to store workplace sizes per workplace per worker group
    # & workplace information (covid secure, workplace open)
    workplace_sizes = Array{Array{Int64,1},1}(undef, workertypes)
    workplace_info = Array{Array{workplace_params,1},1}(undef, workertypes)

    # Initialise ID vars
    workplace_ids = Array{Array{Int64,1},1}(undef, workertypes)  # representing workplace id per person per group
    workertype_ids = Int64[]  # representing workertype id per person

    # Initialise array to store ids of workers within each workplace
    nodes_by_workplace = Array{Array{Array{Int64,1},1},1}(undef,workertypes)

    ## Now generate workplaces within each worker type
    for worker_grp_idx = 1:workertypes
        # initialise arrays for this worker type
        workplace_info[worker_grp_idx] =  workplace_params[]
        workplace_sizes[worker_grp_idx] = Int64[]
        workplace_ids[worker_grp_idx] = Int64[]

        # Initialise array to store ids of workers within each workplace
        nodes_by_workplace[worker_grp_idx] = Array[]

        # Set up tracking counters
        workplace_count = 0   # track number of workplaces in this worker type
        worker_count = 0        # track number of workers assigned to workplaces

        # Add workplaces until worker_count exceeds target number of workers in that sector
        while worker_count < worker_numbers[worker_grp_idx]

            # Update workplace number counter
            workplace_count += 1

            # generate workplace of random size and add to array
            push!(workplace_sizes[worker_grp_idx], ceil(Int64,abs(randn(rng)*workplace_size_sd[worker_grp_idx]+workplace_size_mean[worker_grp_idx])))
            append!(workplace_ids[worker_grp_idx], repeat([workplace_count], workplace_sizes[worker_grp_idx][end]))

            # Add empty workplace to nodes_by_workplace array
            push!(nodes_by_workplace[worker_grp_idx], Int64[])

            # Update tracking counter values
            worker_count += workplace_sizes[worker_grp_idx][end]
        end

        # remove extra workers
        workplace_sizes[worker_grp_idx][end] = worker_numbers[worker_grp_idx] - sum(workplace_sizes[worker_grp_idx][1:end-1])
        workplace_ids[worker_grp_idx] = workplace_ids[worker_grp_idx][1:sum(workplace_sizes[worker_grp_idx])]

        # randomly shuffle workers between workplaces
        shuffle!(rng,workplace_ids[worker_grp_idx])

        append!(workertype_ids, repeat([worker_grp_idx],sum(workplace_sizes[worker_grp_idx])))

        # Create workplace parameter type for this worker type
        workplace_info[worker_grp_idx] =  Array{workplace_params,1}(undef,workplace_count)
        for workplace_itr = 1:workplace_count
            workplace_info[worker_grp_idx][workplace_itr] = workplace_params(covid_secure = CS_active_flag)
        end

    end
    # randomly shuffle workertypes
    shuffle!(rng,workertype_ids)

    worker_nodes = Array{worker_params,1}(undef,cmax)    ## (atwork, workertype, workplaceid)

    for ii = 1:cmax
        returned_to_work = Int64(rand(rng)<workpercent[workertype_ids[ii]])   # decide if returning to work
        worker_nodes[ii] = worker_params(returned_to_work = returned_to_work,
                                            sector_ID = workertype_ids[ii],
                                            workplace_ID = workplace_ids[workertype_ids[ii]][1])

        # All workers added to workplace, regardless of returned_to_work
        push!(nodes_by_workplace[workertype_ids[ii]][workplace_ids[workertype_ids[ii]][1]], ii)
        deleteat!(workplace_ids[workertype_ids[ii]],1) # remove worker from top of list
    end

    return worker_nodes::Array{worker_params,1},
            workplace_sizes::Array{Array{Int64,1},1},
            workplace_info::Array{Array{workplace_params,1},1},
            nodes_by_workplace::Array{Array{Array{Int64,1},1},1}
end


function generate_contacts(cmax::Int64,
                            endtime::Int64,
                            network_parameters::network_params,
                            workplace_generation_params::workplace_generation_params,
                            nodes_by_workplace::Array{Array{Array{Int64,1},1},1},
                            RNGseed::Int64,
                            rng::MersenneTwister)

    # Inputs:
    # cmax - Number of workers in the system
    # endtime - Number of timesteps per simulation
    # network_parameters
    #   worker nodes - array containing worker info
    #   prob_workertype_contact - probability of making contact with others in worktype (DIFFERENT WORKPLACE)
    #   prob_anyworker_contact - probability of making contact with others in DIFFERENT WORKTYPE
    #   prob_social_contact - probability of making a contact socially
    #   dd_within_workplace - mean degree of workers within workplace
    #   household_size_distribution - distribution of household sizes
    #   RNGseed:Int64 - Value to seed the random number generator

# Outputs:
    # work_contacts_same_workplace, work_contacts_other_workplace, household_contacts, social_contacts
    #           - Vector of vectors with IDs of contacts
    # work_contacts_same_workplace_per_node, work_contacts_other_workplace_per_node,
    #       household_contacts_per_node, social_contacts_per_node - Total number of regular contacts within each setting
    # n_households - Total number of households in the system
    # Outputs with _CS appended. Corresponds to contacts made when workplace is Covid-secure.


    @unpack worker_nodes, workplace_sizes,prob_workertype_contact,prob_anyworker_contact,
        prob_social_contact,dd_within_workplace, household_size_distribution,
        workplace_info, CS_active_flag, workplace_degree_distribution,
        between_workplace_contact_probs, social_group_size_distribution,
        friend_of_friend_prob, max_contacts_social, network_generation_method,
        CS_team_size = network_parameters

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Initialise vector of vectors storing IDs of contacts for each node
    # at workplace (non-covid secure setting)
    work_contacts_same_workplace = Array{Array{Int64,1},1}(undef,cmax)
    work_contacts_other_workplace = Array{Array{Int64,1},1}(undef,cmax)

    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    social_contacts = Array{Array{Int64,1},1}(undef,cmax)
    household_contacts = Array{Array{Int64,1},1}(undef,cmax)

    # If CS settings are active,
    # initialise vector of vectors storing IDs of contacts for each node
    # at workplace (covid-secure setting) & vector giving total contacts made by node
    if CS_active_flag == true
        work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
        work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
    end

    for ii = 1:cmax
        work_contacts_same_workplace[ii] = Int64[]
        work_contacts_other_workplace[ii] = Int64[]
        social_contacts[ii] = Int64[]
        household_contacts[ii] = Int64[]

        if CS_active_flag == true
            work_contacts_same_workplace_CS[ii] = Int64[]
        end
    end

    # Initialise vectors giving total contacts made by node
    work_contacts_same_workplace_per_node = zeros(Int64,cmax)
    work_contacts_other_workplace_per_node = zeros(Int64,cmax)
    social_contacts_per_node = zeros(Int64,cmax)
    household_contacts_per_node = zeros(Int64,cmax)

    # Initialise array to store social groups
    nodes_by_social_group = Array[]

    if network_generation_method == "configuration"

        if CS_active_flag == false

            # Construct links for workplace

            if length(workplace_degree_distribution)==1
                workplace_degree_distribution = repeat(workplace_degree_distribution, workplace_generation_params.workertypes)
            end
            if length(between_workplace_contact_probs)==1
                between_workplace_contact_probs = repeat(between_workplace_contact_probs, workplace_generation_params.workertypes)
            end

            worker_count = 0

            # Cycle through sectors and workplaces
            for worker_grp_idx = 1:workplace_generation_params.workertypes

                # nodes_outside_cluster = collect(Iterators.flatten(nodes_by_workplace[worker_grp_idx]))

                # find nodes in the given sector (these are nodes outside the cluster of a given workplace)
                length_array = 0
                for ii=1:length(nodes_by_workplace[worker_grp_idx])
                    length_array += length(nodes_by_workplace[worker_grp_idx][ii])
                end

                nodes_outside_cluster = zeros(Int64,length_array)
                init_pos = 1
                for ii=1:length(nodes_by_workplace[worker_grp_idx])
                    end_pos = init_pos+length(nodes_by_workplace[worker_grp_idx][ii])-1
                    nodes_outside_cluster[init_pos:end_pos] = nodes_by_workplace[worker_grp_idx][ii]
                    init_pos = end_pos+1
                end

                # set the degree distribution and external contact probablity for that sector
                degree_distribution = workplace_degree_distribution[worker_grp_idx]
                external_contact_prob = between_workplace_contact_probs[worker_grp_idx]

                # iterate through the workplaces in that sector
                for workplace_idx = 1:length(workplace_sizes[worker_grp_idx])

                    # nodes in that workplace are in a cluster
                    nodes_within_cluster = nodes_by_workplace[worker_grp_idx][workplace_idx]

                    if length(nodes_within_cluster) > 0
                        worker_count += length(nodes_within_cluster)
                        configuration_model!(worker_nodes,
                                            external_contact_prob,
                                            degree_distribution,
                                            nodes_within_cluster,
                                            nodes_outside_cluster,
                                            work_contacts_same_workplace,
                                            work_contacts_other_workplace,
                                            work_contacts_same_workplace_per_node,
                                            work_contacts_other_workplace_per_node,
                                            network_parameters,
                                            rng)
                    end
                end
            end
            total_ext = sum(work_contacts_other_workplace_per_node)/2
            total_int = sum(work_contacts_same_workplace_per_node)/2
            println("External workplace contacts: $total_ext")
            println("Internal workplace contacts: $total_int")

        ## CS workplace generation
        else
            CS_workplace_generation!(worker_nodes,
                                    nodes_by_workplace,
                                    work_contacts_same_workplace_CS,
                                    work_contacts_same_workplace_per_node_CS,
                                    network_parameters,
                                    rng)

            total_int = sum(work_contacts_same_workplace_per_node_CS)/2
            println("Internal CS workplace contacts: $total_int")
        end

        degree_distribution = social_group_size_distribution

        configuration_model!(cmax,
                                friend_of_friend_prob,
                                degree_distribution,
                                social_contacts,
                                social_contacts_per_node,
                                max_contacts_social,
                                network_parameters,
                                rng)


    elseif network_generation_method == "ER"

        ER_model!(worker_nodes, cmax, workplace_sizes,
                            dd_within_workplace, prob_workertype_contact, prob_anyworker_contact,
                            work_contacts_same_workplace::Array{Array{Int64,1},1},
                            work_contacts_other_workplace::Array{Array{Int64,1},1},
                            work_contacts_same_workplace_per_node::Array{Int64,1},
                            work_contacts_other_workplace_per_node::Array{Int64,1},
                            rng::MersenneTwister)

        # Construct network links for social
        for ii = 1:(cmax-1)
            for jj = (ii+1):cmax
                ### Social contacts with anyone ###
                if rand(rng) < prob_social_contact
                    # Assign IDs of contact to tuple vectors
                    push!(social_contacts[ii],jj)
                    push!(social_contacts[jj],ii)
                    # Increment number of contacts for nodes ii & jj
                    social_contacts_per_node[ii] += 1
                    social_contacts_per_node[jj] += 1
                end
            end
        end
    end

    # Create household contacts & assign a household ID
    csum_household_size_distribution = cumsum(household_size_distribution) # Get cumulative sum, used in each repeat of while loop
    worker_id = 1     # Initialise variables to be incremented in while loop
    household_id = 1
    while worker_id <= cmax
        household_size = findfirst(x->x>rand(rng), csum_household_size_distribution)  # generate random household size

        # check we don't exceed remaining workers left
        household_size = min(household_size, cmax-worker_id+1)

        # create fully-connected household
        for ii = 0:(household_size-1)
            for jj = 0:(household_size-1)
                # Get household contacts
                if ii != jj
                    push!(household_contacts[worker_id+ii], worker_id+jj)
                    household_contacts_per_node[worker_id+ii] += 1
                end
            end

            # Assign household ID to node ii
            worker_nodes[worker_id+ii].household_ID = household_id
        end

        worker_id += household_size
        household_id += 1
    end

    # Assign number of households in total to output variable
    n_households = household_id

    # Define what is returned from the function
    if CS_active_flag == true
        contacts = contacts_struct(cmax=cmax,endtime=endtime,
        work_contacts_same_workplace=work_contacts_same_workplace,
        work_contacts_other_workplace=work_contacts_other_workplace,
        social_contacts_per_node = social_contacts_per_node,
        work_contacts_same_workplace_CS = work_contacts_same_workplace_CS,
        work_contacts_same_workplace_per_node_CS = work_contacts_same_workplace_per_node_CS)
    else
        contacts = contacts_struct(cmax=cmax,endtime=endtime,work_contacts_same_workplace=work_contacts_same_workplace,
        work_contacts_other_workplace=work_contacts_other_workplace,
        social_contacts_per_node = social_contacts_per_node)
    end

    return contacts::contacts_struct,
        social_contacts::Array{Array{Int64,1},1},
        household_contacts::Array{Array{Int64,1},1},
        work_contacts_same_workplace_per_node::Array{Int64,1},
        work_contacts_other_workplace_per_node::Array{Int64,1},
        household_contacts_per_node::Array{Int64,1},
        n_households::Int64
end


#### Covid-secure workplace generation ####
function CS_workplace_generation!(worker_nodes::Array{worker_params,1},
                                    nodes_by_workplace::Array{Array{Array{Int64,1},1},1},
                                    work_contacts_same_workplace_CS::Array{Array{Int64,1},1},
                                    work_contacts_same_workplace_per_node_CS::Array{Int64,1},
                                    network_parameters::network_params,
                                    rng::MersenneTwister)

    @unpack workplace_sizes, CS_team_size = network_parameters

    workertypes = length(nodes_by_workplace)

    # Cycle through sectors and workplaces
    for worker_grp_idx = 1:workertypes

        n_workplaces = length(nodes_by_workplace[worker_grp_idx])

        # iterate through the workplaces in that sector
        for workplace_idx = 1:n_workplaces

            # Find nodes in this workplace
            nodes_within_workplace = nodes_by_workplace[worker_grp_idx][workplace_idx]
            # Number of nodes in workplace
            n_nodes_within_workplace = workplace_sizes[worker_grp_idx][workplace_idx]
            # Number of isolated, fully connected groups of size CS_team_size
            n_CS_groups = floor(Int64, n_nodes_within_workplace/CS_team_size)

            # Create fully connected groups
            for CS_group_idx = 1:n_CS_groups
                current_group = nodes_within_workplace[((CS_group_idx-1)*CS_team_size+1):(CS_group_idx*CS_team_size)]
                for worker_idx = 1:CS_team_size
                    node_id = current_group[worker_idx]
                    work_contacts_same_workplace_CS[node_id] = current_group[1:end .!= worker_idx]
                    work_contacts_same_workplace_per_node_CS[node_id] = CS_team_size - 1
                end
            end

            # Group together any remaining workers in workplace
            workers_remaining = Int64(n_nodes_within_workplace%CS_team_size)

            if workers_remaining > 1
                current_group = nodes_within_workplace[(n_CS_groups*CS_team_size+1):end]
                for worker_idx = 1:workers_remaining
                    node_id = current_group[worker_idx]
                    work_contacts_same_workplace_CS[node_id] = current_group[1:end .!= worker_idx]
                    work_contacts_same_workplace_per_node_CS[node_id] = workers_remaining - 1
                end
            end
        end
    end
    return nothing
end



#### Workplace version ####
function configuration_model!(worker_nodes::Array{worker_params,1},
                                external_contact_prob::Float64,
                                degree_distribution::Distribution,
                                nodes_within_cluster::Array{Int64,1},
                                nodes_outside_cluster::Array{Int64,1},
                                work_contacts_same_workplace::Array{Array{Int64,1},1},
                                work_contacts_other_workplace::Array{Array{Int64,1},1},
                                work_contacts_same_workplace_per_node::Array{Int64,1},
                                work_contacts_other_workplace_per_node::Array{Int64,1},
                                network_parameters::network_params,
                                rng::MersenneTwister)


    n_nodes = length(nodes_within_cluster)
    n_nodes_external = length(nodes_outside_cluster)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Decrease by proportion of contacts made external to workplace
    # Limit maximum number of contacts to workplace_size - 1
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, (edges_per_node[node_id])/(1+external_contact_prob))
        edges_per_node[node_id] = min((n_nodes - 1), edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    n_stubs = half_edges[end]

    edges_within_group = zeros(Int64, n_nodes, n_nodes)
    [edges_within_group[i,i] = 1 for i=1:n_nodes]

    n_nodes_external_remaining = n_nodes_external

    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # Contact made external to cluster
        if (rand(rng) < external_contact_prob) & (n_nodes_external_remaining > n_nodes)

            node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)

            # Don't allow links within cluster or repeated links
            while (nodes_outside_cluster[node_id2] in work_contacts_other_workplace[nodes_within_cluster[node_id1]]) || (node_id2 in nodes_within_cluster)
                node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)
            end

            if (worker_nodes[nodes_within_cluster[node_id1]].returned_to_work == 1) &&
                (worker_nodes[nodes_outside_cluster[node_id2]].returned_to_work == 1)
                push!(work_contacts_other_workplace[nodes_within_cluster[node_id1]], nodes_outside_cluster[node_id2])
                push!(work_contacts_other_workplace[nodes_outside_cluster[node_id2]], nodes_within_cluster[node_id1])

                work_contacts_other_workplace_per_node[nodes_within_cluster[node_id1]] += 1
                work_contacts_other_workplace_per_node[nodes_outside_cluster[node_id2]] += 1
            end

            n_nodes_external_remaining -= 1

            for ii=node_id1:n_nodes
                half_edges[ii] -= 1
            end

        # Contact made within cluster
        else
            # This ensures no self-links or repeated edges within cluster
            nodes_remaining = findall(edges_within_group[node_id1,:] .== 0)
            half_edges_remaining = cumsum(edges_per_node[nodes_remaining])
            # removals = findall(half_edges_remaining.==0)
            # deleteat!(nodes_remaining,removals)
            # deleteat!(half_edges_remaining,removals)
            nodes_remaining = nodes_remaining[half_edges_remaining.>0]
            half_edges_remaining = half_edges_remaining[half_edges_remaining.>0]

            # Contact not made (half edge lost)
            if length(nodes_remaining) < 1

                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end
            else
                stub_id2 = round(Int64, rand(rng)*(half_edges_remaining[end]-1) + 1)
                node_id2 = nodes_remaining[findfirst(x -> x >= stub_id2, half_edges_remaining)]

                if (worker_nodes[nodes_within_cluster[node_id1]].returned_to_work == 1) &&
                    (worker_nodes[nodes_within_cluster[node_id2]].returned_to_work == 1)
                    push!(work_contacts_same_workplace[nodes_within_cluster[node_id1]], nodes_within_cluster[node_id2])
                    push!(work_contacts_same_workplace[nodes_within_cluster[node_id2]], nodes_within_cluster[node_id1])

                    work_contacts_same_workplace_per_node[nodes_within_cluster[node_id1]] += 1
                    work_contacts_same_workplace_per_node[nodes_within_cluster[node_id2]] += 1
                end

                edges_within_group[node_id1,node_id2] += 1

                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end
                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end
            end
        end
    end

end

#### Social version ####
function configuration_model!(external_contact_prob::Float64,
                                degree_distribution::Distribution,
                                nodes_within_cluster::Array{Int64,1},
                                nodes_outside_cluster::Array{Int64,1},
                                social_contacts::Array{Array{Int64,1},1},
                                social_contacts_per_node::Array{Int64,1},
                                max_contacts_social::Int64,
                                rng::MersenneTwister)

    n_nodes = length(nodes_within_cluster)
    n_nodes_external = length(nodes_outside_cluster)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Check if each node has already made contacts from other clusters
    # Limit maximum number of contacts to max_contacts_social
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, edges_per_node[node_id])
        edges_per_node[node_id] -= length(social_contacts[nodes_within_cluster[node_id]])
        edges_per_node[node_id] = min(max_contacts_social, edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    n_stubs = half_edges[end]

    edges_within_group = zeros(Int64, n_nodes, n_nodes)
    [edges_within_group[i,i] = 1 for i=1:n_nodes]

    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # This ensures no self-links or repeated edges within cluster
        nodes_remaining = findall(edges_within_group[node_id1,:] .== 0)
        half_edges_remaining = cumsum(edges_per_node[nodes_remaining])
        nodes_remaining = nodes_remaining[half_edges_remaining.>0]
        half_edges_remaining = half_edges_remaining[half_edges_remaining.>0]

        # Contact made external to cluster
        if (length(nodes_remaining) < 1) || (rand(rng) < external_contact_prob)

            node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)

            # Don't allow links within cluster
            while node_id2 in nodes_within_cluster
                node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)
            end

            push!(social_contacts[nodes_within_cluster[node_id1]], nodes_outside_cluster[node_id2])
            push!(social_contacts[nodes_outside_cluster[node_id2]], nodes_within_cluster[node_id1])

            social_contacts_per_node[nodes_within_cluster[node_id1]] += 1
            social_contacts_per_node[nodes_outside_cluster[node_id2]] += 1


            for ii=node_id1:n_nodes
                half_edges[ii] -= 1
            end

        # Contact made within cluster
        else
            stub_id2 = round(Int64, rand(rng)*(half_edges_remaining[end]-1) + 1)
            node_id2 = nodes_remaining[findfirst(x -> x >= stub_id2, half_edges_remaining)]

            push!(social_contacts[nodes_within_cluster[node_id1]], nodes_within_cluster[node_id2])
            push!(social_contacts[nodes_within_cluster[node_id2]], nodes_within_cluster[node_id1])

            social_contacts_per_node[nodes_within_cluster[node_id1]] += 1
            social_contacts_per_node[nodes_within_cluster[node_id2]] += 1

            edges_within_group[node_id1,node_id2] += 1


            for ii=node_id1:n_nodes
                half_edges[ii] -= 1
            end
            for ii=node_id2:n_nodes
                half_edges[ii] -= 1
            end
        end
    end

end


#### Social version (friends-of-friends model) ####
function configuration_model!(n_nodes::Int64,
                                friend_of_friend_prob::Float64,
                                degree_distribution::Distribution,
                                social_contacts::Array{Array{Int64,1},1},
                                social_contacts_per_node::Array{Int64,1},
                                max_contacts_social::Int64,
                                network_parameters::network_params,
                                rng::MersenneTwister)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Limit maximum number of contacts to max_contacts_social
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, edges_per_node[node_id])
        edges_per_node[node_id] = min(max_contacts_social, edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    if half_edges[end] % 2 != 0
        half_edges[end] += 1
    end

    n_stubs = half_edges[end]

    edges_remaining_per_node = copy(edges_per_node)
    fof = Int64[]
    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # fof = collect(Iterators.flatten(social_contacts[social_contacts[node_id1]]))
        fof = find_fof(social_contacts,social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)

        # Contact made with friend of friend
        if (length(fof) > 0) & (rand(rng) < friend_of_friend_prob)

            node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

        # Contact made with not friend of friend
        else

            stub_id2 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id2 = findfirst(x -> x >= stub_id2, half_edges)

        end

        push!(social_contacts[node_id1], node_id2)
        push!(social_contacts[node_id2], node_id1)

        social_contacts_per_node[node_id1] += 1
        social_contacts_per_node[node_id2] += 1

        edges_remaining_per_node[node_id1] -= 1
        edges_remaining_per_node[node_id2] -= 1


        for ii=node_id1:n_nodes
            half_edges[ii] -= 1
        end
        for ii=node_id2:n_nodes
            half_edges[ii] -= 1
        end

    end

end


function ER_model!(worker_nodes, cmax, workplace_sizes,
                    dd_within_workplace, prob_workertype_contact, prob_anyworker_contact,
                    work_contacts_same_workplace::Array{Array{Int64,1},1},
                    work_contacts_other_workplace::Array{Array{Int64,1},1},
                    work_contacts_same_workplace_per_node::Array{Int64,1},
                    work_contacts_other_workplace_per_node::Array{Int64,1},
                    rng::MersenneTwister)

    for ii = 1:(cmax-1)

        returned_to_work::Int64 = worker_nodes[ii].returned_to_work
        worker_grp_idx::Int64 = worker_nodes[ii].sector_ID
        workplace_idx::Int64 = worker_nodes[ii].workplace_ID

        for jj = (ii+1):cmax

            # If returned to work, WORK contacts consist of contacts within workplace,
            # + lower level contacts with workers of same worker group
            # Else, if WFH, no work contacts

            ### On days when at work, increased contacts with others at work  ###
            if (returned_to_work == 1) & (worker_nodes[jj].returned_to_work==1)  # both returned to work

                ## For workers in the same group and workplace, edges form according to ER graph
                if (worker_nodes[jj].sector_ID == worker_grp_idx) & (worker_nodes[jj].workplace_ID == workplace_idx)
                    if rand(rng) < dd_within_workplace[worker_grp_idx]/(workplace_sizes[worker_grp_idx][workplace_idx] - 1) ## ER component
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_same_workplace[ii],jj)
                        push!(work_contacts_same_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_same_workplace_per_node[ii] += 1
                        work_contacts_same_workplace_per_node[jj] += 1
                    end

                ## same worker group, different workplace
                elseif (worker_nodes[jj].sector_ID == worker_grp_idx)

                    if rand(rng) < prob_workertype_contact[worker_grp_idx]
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_other_workplace[ii],jj)
                        push!(work_contacts_other_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_other_workplace_per_node[ii] += 1
                        work_contacts_other_workplace_per_node[jj] += 1
                    end

                ## different worker group
                else
                    if rand(rng) < prob_anyworker_contact
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_other_workplace[ii],jj)
                        push!(work_contacts_other_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_other_workplace_per_node[ii] += 1
                        work_contacts_other_workplace_per_node[jj] += 1
                    end
                end
            end
        end
    end

end



# Premake dynamic worker contacts,
# to be loaded in ahead of simulation
function generate_dynamic_worker_contacts(RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            worker_nodes::Array{worker_params,1},
                                            dynamic_work_dd::Array{Distribution,1},
                                            max_contacts_work_dynamic::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# worker nodes - array with entry per worker
# dynamic_conts_mean, dynamic_conts_sd - Distribution properties for dynamic worker contacts

# Outputs:
# dynamic_worker_contacts - Per node, a record of dynamic worker contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_worker_contacts = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for node_itr = 1:cmax
        if worker_nodes[node_itr].returned_to_work==1 # Add dynamic links, if returned to work
            # Get dynamic worker group type for node_itr
            node_dynamic_grp_ID = worker_nodes[node_itr].sector_ID

            for time_itr = 1:endtime
                # Generate number of dynamic contacts from appropriate distribution
                gg = round(Int64, Distributions.rand(rng, dynamic_work_dd[node_dynamic_grp_ID]))
                # Limit to specified maximum
                gg = min(gg, max_contacts_work_dynamic)
                # If dynamic worker contacts made, assign to output variable
                if gg > 0
                    dynamic_worker_contacts[time_itr,node_itr] = zeros(gg)

                    # Generate required number of contacts
                    for contact_itr = 1:gg
                        gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        while gg1 == node_itr # Redraw if returned the index node
                            gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        end
                        dynamic_worker_contacts[time_itr,node_itr][contact_itr] = gg1
                    end
                else # No dynamic worker contacts made on given timestep
                    dynamic_worker_contacts[time_itr,node_itr] = Int64[]
                end
            end
        end
    end

    return dynamic_worker_contacts::Array{Array{Int64,1},2}

end

function find_fof(contacts::Array{Array{Int64,1},1},social_contacts_of_node::Array{Int64,1},
    edges_remaining_per_node::Array{Float64,1},node_id1::Int64,fof::Array{Int64,1})

    # # Find groups of friends already meeting today with edges remaining
    # length_array = 0
    # for ii=1:length(contacts[social_contacts_of_node])
    #     length_array += length(contacts[social_contacts_of_node][ii])
    # end
    #
    # fof = zeros(Int64,length_array)
    # init_pos = 1
    # for ii=1:length(contacts[social_contacts_of_node])
    #     end_pos = init_pos+length(contacts[social_contacts_of_node][ii])-1
    #     fof[init_pos:end_pos] = contacts[social_contacts_of_node][ii]
    #     init_pos = end_pos+1
    # end

    fof = collect(Iterators.flatten(contacts[social_contacts_of_node]))

    # # Remove self and anyone already friends with
    # fof = fof[edges_remaining_per_node[fof].>0]
    # fof = fof[fof .!= node_id1]
    # fof = fof[[fof[i] ∉ social_contacts_of_node for i=1:length(fof)]]

    # Remove self and anyone already friends with and those without edges remaining
    # This version doesn't keep duplicates, so is no more likely to make an edge
    # with a node that is the friend of multiple friends
    # makes fewer allocations and is a lot faster
    setdiff!(fof,social_contacts_of_node)
    setdiff!(fof,[node_id1])
    removals = findall(edges_remaining_per_node[fof].==0)
    deleteat!(fof,removals)

    return fof
end


# Premake daily social contacts,
# to be loaded in ahead of simulation
# CONFIGURATION MODEL VERSION
function generate_social_contacts_each_day(rng::MersenneTwister,
                                    RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_contacts_per_node::Array{Int64,1},
                                            social_workday_dd::Distribution,
                                            social_nonworkday_dd::Distribution,
                                            cluster_social_contacts::Bool)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# social_contacts_per_node - Total amount of social contacts each individual has (entry per individual)
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Configuration model
    """
    n_nodes = cmax

    edges_per_node_workday = Distributions.rand(rng, social_workday_dd, n_nodes)
    edges_per_node_nonworkday = Distributions.rand(rng, social_workday_dd, n_nodes)

    # Round degree to nearest whole number
    # Limit maximum number of contacts to size of friendship group
    for node_id = 1:n_nodes
        friendship_group_size = social_contacts_per_node[node_id]
        edges_per_node_workday[node_id] = round(Int64, edges_per_node_workday[node_id])
        edges_per_node_workday[node_id] = min(friendship_group_size, edges_per_node_workday[node_id])
        edges_per_node_nonworkday[node_id] = round(Int64, edges_per_node_nonworkday[node_id])
        edges_per_node_nonworkday[node_id] = min(friendship_group_size, edges_per_node_nonworkday[node_id])
        # Initialise array to store daily contacts
        for time_itr = 1:endtime
            workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            nonworkday_social_contacts_by_day[time_itr,node_id] = Int64[]
        end
    end

    fof = Int64[]
    for time_itr = 1:endtime
        println(time_itr)

        """
        WORKDAY
        """
        half_edges = cumsum(edges_per_node_workday)

        edges_remaining_per_node = copy(edges_per_node_workday)
        while half_edges[end] > 1

            stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id1 = findfirst(x -> x >= stub_id1, half_edges)

            if cluster_social_contacts == true
                # Find groups of friends already meeting today with edges remaining
                # fof = collect(Iterators.flatten(workday_social_contacts_by_day[time_itr,social_contacts[node_id1]]))
                fof = find_fof(workday_social_contacts_by_day,social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)
            else
                fof = Int64[]
            end

            # Contact made with group of friends already meeting
            if length(fof) > 0

                node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

                push!(workday_social_contacts_by_day[time_itr,node_id1], node_id2)
                push!(workday_social_contacts_by_day[time_itr,node_id2], node_id1)

                edges_remaining_per_node[node_id1] -= 1
                edges_remaining_per_node[node_id2] -= 1


                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end

                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end

            # Contact made with anyone from friend group
            else
                friend_group = social_contacts[node_id1]

                # Remove anyone already meeting today
                friend_group = friend_group[[friend_group[i] ∉ workday_social_contacts_by_day[time_itr,node_id1] for i=1:length(friend_group)]]
                # Remove those with no half edges remaining
                friend_group = friend_group[edges_remaining_per_node[friend_group].>0]

                # If possible, connect to random friend
                if length(friend_group) > 0
                    node_id2 = friend_group[ceil(Int64, rand(rng)*length(friend_group))]

                    push!(workday_social_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(workday_social_contacts_by_day[time_itr,node_id2], node_id1)

                    edges_remaining_per_node[node_id1] -= 1
                    edges_remaining_per_node[node_id2] -= 1


                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end

                    for ii=node_id2:n_nodes
                        half_edges[ii] -= 1
                    end

                # Otherwise discard half-edge
                else
                    edges_remaining_per_node[node_id1] -= 1

                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end
                end
            end
        end

        """
        NON-WORKDAY
        """
        half_edges = cumsum(edges_per_node_nonworkday)

        edges_remaining_per_node = copy(edges_per_node_nonworkday)
        while half_edges[end] > 1

            stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id1 = findfirst(x -> x >= stub_id1, half_edges)

            # fof = collect(Iterators.flatten(nonworkday_social_contacts_by_day[time_itr,social_contacts[node_id1]]))
            fof = find_fof(view(nonworkday_social_contacts_by_day[time_itr,:]),social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)

            # Contact made with group of friends already meeting
            if (length(fof) > 0) & (cluster_social_contacts == true)

                node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

                push!(nonworkday_social_contacts_by_day[time_itr,node_id1], node_id2)
                push!(nonworkday_social_contacts_by_day[time_itr,node_id2], node_id1)

                edges_remaining_per_node[node_id1] -= 1
                edges_remaining_per_node[node_id2] -= 1


                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end

                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end

            # Contact made with anyone from friend group
            else
                friend_group = social_contacts[node_id1]

                # Remove anyone already meeting today
                friend_group = friend_group[[friend_group[i] ∉ nonworkday_social_contacts_by_day[time_itr,node_id1] for i=1:length(friend_group)]]
                # Remove those with no half edges remaining
                friend_group = friend_group[edges_remaining_per_node[friend_group].>0]

                # If possible, connect to random friend
                if length(friend_group) > 0
                    node_id2 = friend_group[ceil(Int64, rand(rng)*length(friend_group))]

                    push!(nonworkday_social_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(nonworkday_social_contacts_by_day[time_itr,node_id2], node_id1)

                    edges_remaining_per_node[node_id1] -= 1
                    edges_remaining_per_node[node_id2] -= 1


                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end

                    for ii=node_id2:n_nodes
                        half_edges[ii] -= 1
                    end

                # Otherwise discard half-edge
                else
                    edges_remaining_per_node[node_id1] -= 1

                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end
                end
            end
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end

### Cluster version ###
function generate_social_contacts_each_day(rng::MersenneTwister,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_workday_dd::Distribution,
                                            social_nonworkday_dd::Distribution)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    for time_itr = 1:endtime
        # println(time_itr)
        """
        WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) > 0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_workday_dd)

            # Round and limit to number of unassigned friends
            cluster_size = round(Int64, min(cluster_size, length(friend_list)))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                workday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    workday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    if jj>length(other_friends)
                        workday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end

                    # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # spots_remaining = cluster_size - (ii-1)
                    #
                    # # friends_of_friends_in_cluster = friends_in_cluster[[friends_in_cluster[i] ∈ social_contacts[friend_id] for i=1:cluster_size]]
                    # # # Assign intersection to today's contacts
                    # # # if length(friends_of_friends_in_cluster) > 0
                    # #     workday_social_contacts_by_day[time_itr,friend_id] = friends_of_friends_in_cluster
                    # # # else
                    # # #     workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    # # # end
                    # # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # # spots_remaining = cluster_size - length(friends_of_friends_in_cluster)
                    # # Assign remaining spots to other friends at random
                    # other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    # spots_remaining = min(spots_remaining, length(other_friends))
                    # shuffle!(rng,other_friends)
                    # for spot_itr = 1:spots_remaining
                    #     push!(workday_social_contacts_by_day[time_itr,friend_id], other_friends[spot_itr])
                    # end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end

        """
        NON-WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) >0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_nonworkday_dd)

            # Round and limit to number of unassigned friends
            cluster_size = round(Int64, min(cluster_size, length(friend_list)))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                nonworkday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    nonworkday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            nonworkday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            nonworkday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    # if there weren't enough other friends to fill the contacts, remove extra zeros
                    if jj>length(other_friends)
                        nonworkday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end

                    # # Find intersection of friend group with current cluster
                    # friends_of_friends_in_cluster = friends_in_cluster[[friends_in_cluster[i] ∈ social_contacts[friend_id] for i=1:cluster_size]]
                    # # Assign intersection to today's contacts
                    # if length(friends_of_friends_in_cluster) > 0
                    #     nonworkday_social_contacts_by_day[time_itr,friend_id] = friends_of_friends_in_cluster
                    # else
                    #     nonworkday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    # end
                    # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # spots_remaining = cluster_size - length(friends_of_friends_in_cluster)
                    # # Assign remaining spots to other friends at random
                    # other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    # spots_remaining = min(spots_remaining, length(other_friends))
                    # shuffle!(other_friends)
                    # for spot_itr = 1:spots_remaining
                    #     push!(nonworkday_social_contacts_by_day[time_itr,friend_id], other_friends[spot_itr])
                    # end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                nonworkday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


#### Version to assess role of static vs dynamic social contacts
function generate_social_contacts_each_day(rng::MersenneTwister,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_workday_dd::Distribution,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = 1
    while time_itr <= endtime
        # println(time_itr)
        """
        WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) > 0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_workday_dd)

            # Round and limit to number of unassigned friends or group limit
            cluster_size = round(Int64, min(cluster_size, length(friend_list), group_limit))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                workday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    workday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    if jj>length(other_friends)
                        workday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


#### Toy model for rule of 6 analysis
function generate_social_contacts_each_day(rng::MersenneTwister,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            n_groups_per_day_distribution::Distribution,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = 1
    while time_itr <= endtime

        n_groups_per_node = zeros(Int64, cmax)
        for node_itr = 1:cmax
            n_groups_per_node[node_itr] = Distributions.rand(rng, n_groups_per_day_distribution)
            workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
        end

        while sum(n_groups_per_node) > 0

            nodes_remaining = findall(n_groups_per_node.>0)

            current_group = shuffle(nodes_remaining)[1:min(group_limit,length(nodes_remaining))]

            group_size = length(current_group)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end

            n_groups_per_node[current_group] .-= 1

        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


function generate_social_contacts_each_day(rng::MersenneTwister,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            contacts_per_day::Int64,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = 1
    while time_itr <= endtime

        n_groups_per_node = ones(Int64, cmax)*round(Int64,contacts_per_day/group_limit)
        for node_itr = 1:cmax
            workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
        end

        while sum(n_groups_per_node) > 0

            nodes_remaining = findall(n_groups_per_node.>0)

            current_group = shuffle(nodes_remaining)[1:min(group_limit,length(nodes_remaining))]

            group_size = length(current_group)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end

            n_groups_per_node[current_group] .-= 1

        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end



# # Premake dynamic worker contacts,
# # to be loaded in ahead of simulation
# function generate_social_contacts_each_day(rng::MersenneTwister,
#                                             RNGseed::Int64,
#                                             cmax::Int64,
#                                             endtime::Int64,
#                                             social_contacts::Array{Array{Int64,1},1},
#                                             social_contacts_per_node::Array{Int64,1},
#                                             social_workday_dd::Distribution,
#                                             social_nonworkday_dd::Distribution)
# # Inputs:
# # RNGseed - Seed the random number generator
# # cmax - Number of nodes in the system
# # endtime - Number of timesteps simulation will run for
# # social_contacts - array of arrays containing all possible social contacts for each individual
# # social_contacts_per_node - Total amount of social contacts each individual has (entry per individual)
# # n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts
#
# # Outputs:
# # workday_social_contacts_by_day, nonworkday_social_contacts_by_day
# #       - Per node, a record of social contacts made on each day
#
#     """
#     Set the RNG
#     """
#     rng = MersenneTwister(RNGseed)
#
#     """
#     Initialise vector of vectors storing IDs of contacts for each node
#     """
#     workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
#     nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
#
#     """
#     Iterate over all nodes
#     For those potentially have social contacts,
#     assign a sample of those contacts for each timestep
#     """
#     for node_itr = 1:cmax
#
#         # Add in social contacts if possible
#         if social_contacts_per_node[node_itr] > 0
#
#             # Social contacts may differ each day
#             for time_itr = 1:endtime
#
#                 # Draw random number of social links to be made today, limited to whole friend group
#                 n_social_workday = min(round(Int64, Distributions.rand(rng, social_workday_dd)), social_contacts_per_node[node_itr])
#                 n_social_nonworkday = min(round(Int64, Distributions.rand(rng, social_nonworkday_dd)), social_contacts_per_node[node_itr])
#
#                 if n_social_workday >0
#                     workday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_workday)
#
#                     for wd_social_it = 1:n_social_workday
#                         wd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
#                         workday_social_contacts_by_day[time_itr,node_itr][wd_social_it] = social_contacts[node_itr][wd_contact_idx]
#                     end
#                 else
#                     workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
#                 end
#
#                 if n_social_nonworkday >0
#                     nonworkday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_nonworkday)
#
#                     for nwd_social_it = 1:n_social_nonworkday
#                         nwd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
#                         nonworkday_social_contacts_by_day[time_itr,node_itr][nwd_social_it] = social_contacts[node_itr][nwd_contact_idx]
#                     end
#                 else
#                     nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
#                 end
#             end
#
#         end
#     end
#
#     return workday_social_contacts_by_day::Array{Array{Int64,1},2},
#             nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
# end


"""
Functions for use with Erdos-Renyi network construction
"""
# ASSUMPTIONS:
# Contacts are split into 3 groups: work, social and household
# Work contacts occur within workplaces (common), between same workertypes (less common) and between any other worker (uncommon) on every work day
# If working from home, there are NO work contacts
# Social contacts can be made with anyone, but are NOT static - instead take a subset each day (can be a larger subset on weekends)
# Household contacts are within a house and occur everyday
# During '3 days off', workers are NOT working (i.e. more sociable)
function generate_contacts_ER_within_workplace(cmax::Int64,endtime::Int64,network_parameters::network_params,RNGseed::Int64)

   # Inputs:
   # cmax - Number of workers in the system
   # endtime - Number of timesteps per simulation
   # network_parameters
   #   worker nodes - array containing worker info
   #   prob_workertype_contact - probability of making contact with others in worktype (DIFFERENT WORKPLACE)
   #   prob_anyworker_contact - probability of making contact with others in DIFFERENT WORKTYPE
   #   prob_social_contact - probability of making a contact socially
   #   dd_within_workplace - mean degree of workers within workplace
   #   household_size_distribution - distribution of household sizes
   #   RNGseed:Int64 - Value to seed the random number generator

# Outputs:
   # work_contacts_same_workplace, work_contacts_other_workplace, household_contacts, social_contacts
   #           - Vector of vectors with IDs of contacts
   # work_contacts_same_workplace_per_node, work_contacts_other_workplace_per_node,
   #       household_contacts_per_node, social_contacts_per_node - Total number of regular contacts within each setting
   # n_households - Total number of households in the system
   # Outputs with _CS appended. Corresponds to contacts made when workplace is Covid-secure.


   @unpack worker_nodes, workplace_sizes,prob_workertype_contact,prob_anyworker_contact,
       prob_social_contact,dd_within_workplace, household_size_distribution,
       workplace_info, CS_active_flag = network_parameters

   # Set the RNG
   rng = MersenneTwister(RNGseed)

   # Initialise vector of vectors storing IDs of contacts for each node
   # at workplace (non-covid secure setting)
   work_contacts_same_workplace = Array{Array{Int64,1},1}(undef,cmax)
   work_contacts_other_workplace = Array{Array{Int64,1},1}(undef,cmax)

   # Initialise vector of vectors storing IDs of contacts for each node in social
   # and household settings
   social_contacts = Array{Array{Int64,1},1}(undef,cmax)
   household_contacts = Array{Array{Int64,1},1}(undef,cmax)

   # If CS settings are active,
   # initialise vector of vectors storing IDs of contacts for each node
   # at workplace (covid-secure setting) & vector giving total contacts made by node
   if CS_active_flag == true
       work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
       work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
   end

   for ii = 1:cmax
       work_contacts_same_workplace[ii] = Int64[]
       work_contacts_other_workplace[ii] = Int64[]
       social_contacts[ii] = Int64[]
       household_contacts[ii] = Int64[]

       if CS_active_flag == true
           work_contacts_same_workplace_CS[ii] = Int64[]
       end
   end

   # Initialise vectors giving total contacts made by node
   work_contacts_same_workplace_per_node = zeros(Int64,cmax)
   work_contacts_other_workplace_per_node = zeros(Int64,cmax)
   social_contacts_per_node = zeros(Int64,cmax)
   household_contacts_per_node = zeros(Int64,cmax)

   # Construct network links for social, household, & workplace contacts
   # Workplace contacts done for non-covid secure setting
   for ii = 1:(cmax-1)

       returned_to_work::Int64 = worker_nodes[ii].returned_to_work
       worker_grp_idx::Int64 = worker_nodes[ii].sector_ID
       workplace_idx::Int64 = worker_nodes[ii].workplace_ID

       for jj = (ii+1):cmax

           # SOCIAL contacts consist of lower level contacts with anyone
           # If returned to work, WORK contacts consist of contacts within workplace,
           # + lower level contacts with workers of same worker group
           # Else, if WFH, no work contacts

           ### Social contacts with anyone ###
           if rand(rng) < prob_social_contact
               # Assign IDs of contact to tuple vectors
               push!(social_contacts[ii],jj)
               push!(social_contacts[jj],ii)
               # Increment number of contacts for nodes ii & jj
               social_contacts_per_node[ii] += 1
               social_contacts_per_node[jj] += 1
           end

           ### On days when at work, increased contacts with others at work  ###
           if (returned_to_work == 1) & (worker_nodes[jj].returned_to_work==1)  # both returned to work

               ## For workers in the same group and workplace, edges form according to ER graph
               if (worker_nodes[jj].sector_ID == worker_grp_idx) & (worker_nodes[jj].workplace_ID == workplace_idx)
                   if rand(rng) < dd_within_workplace[worker_grp_idx]/(workplace_sizes[worker_grp_idx][workplace_idx] - 1) ## ER component
                       # Assign IDs of contact to tuple vectors
                       push!(work_contacts_same_workplace[ii],jj)
                       push!(work_contacts_same_workplace[jj],ii)

                       # Increment number of contacts for nodes ii & jj
                       work_contacts_same_workplace_per_node[ii] += 1
                       work_contacts_same_workplace_per_node[jj] += 1
                   end

               ## same worker group, different workplace
           elseif (worker_nodes[jj].sector_ID == worker_grp_idx)

                   if rand(rng) < prob_workertype_contact[worker_grp_idx]
                       # Assign IDs of contact to tuple vectors
                       push!(work_contacts_other_workplace[ii],jj)
                       push!(work_contacts_other_workplace[jj],ii)

                       # Increment number of contacts for nodes ii & jj
                       work_contacts_other_workplace_per_node[ii] += 1
                       work_contacts_other_workplace_per_node[jj] += 1
                   end

               ## different worker group
               else
                   if rand(rng) < prob_anyworker_contact
                       # Assign IDs of contact to tuple vectors
                       push!(work_contacts_other_workplace[ii],jj)
                       push!(work_contacts_other_workplace[jj],ii)

                       # Increment number of contacts for nodes ii & jj
                       work_contacts_other_workplace_per_node[ii] += 1
                       work_contacts_other_workplace_per_node[jj] += 1
                   end
               end
           end
       end
   end

   # If covid-secure workplaces in use,
   # construct in-workplace contacts assuming a covid secure setting
   if CS_active_flag == true
       node_CS_assigned = zeros(cmax)
       for node_itr = 1:(cmax-1) # Note, final node will have been accounted for

           # Check if node has already been assigned a workplace "team"
           # If not, go through assignment loop
            if node_CS_assigned[node_itr] == 0

                # Get permitted team size
                relevant_sector = worker_nodes[node_itr].sector_ID
                relevant_workplace = worker_nodes[node_itr].workplace_ID
                permitted_team_size = workplace_info[relevant_sector][relevant_workplace].CS_team_size

                # Load pre-existing workplace contacts
                nonCS_workplace_contacts = work_contacts_same_workplace[node_itr]

                # Check how many workplaces have yet to be allocated
                node_to_be_assigned = convert(Int64,sum(node_CS_assigned[nonCS_workplace_contacts].==0))
                if node_to_be_assigned > 0
                    # Assign workplaces that could be allocated to vector
                    to_be_assigned_workplace_contacts = zeros(Int64,node_to_be_assigned)
                    assignment_idx = 1
                    for assign_itr = 1:length(nonCS_workplace_contacts)
                        ID_to_check = nonCS_workplace_contacts[assign_itr]
                        if node_CS_assigned[ID_to_check] == 0
                            to_be_assigned_workplace_contacts[assignment_idx] = ID_to_check

                            # Increment counter
                            assignment_idx += 1
                        end
                    end
                else # No other workers available to be allocated to this team
                    to_be_assigned_workplace_contacts = Int64[]
                end


                # If number of existing contacts is below (max team size-1),
                # create fully-connected team
                n_to_be_assigned_workplace_contacts = length(to_be_assigned_workplace_contacts)

                # If existing workplace contact amount larger than permitted size,
                # pick a subset, and populate those vectors
                if n_to_be_assigned_workplace_contacts > (permitted_team_size-1)
                    CS_workplace_contact_IDs = sample(rng, to_be_assigned_workplace_contacts, permitted_team_size-1; replace=false, ordered=false)
                else
                    CS_workplace_contact_IDs = copy(to_be_assigned_workplace_contacts)
                end
                n_team_contacts = length(CS_workplace_contact_IDs)

                # Populate workplace_CS contacts for node_itr
                work_contacts_same_workplace_CS[node_itr] = CS_workplace_contact_IDs
                work_contacts_same_workplace_per_node_CS[node_itr] = n_team_contacts


                # Connect each team member to all other team members
                for ii = 1:n_team_contacts

                   # Get ID of team member
                   team_contact_ID = CS_workplace_contact_IDs[ii]

                   # Initialise storage of workplace_CS contacts for node team_contact_ID
                   work_contacts_same_workplace_CS[team_contact_ID] = zeros(Int64,n_team_contacts)

                   # Set number of contacts for node_itr & node team_contact_ID
                   work_contacts_same_workplace_per_node_CS[team_contact_ID] = n_team_contacts

                   # Connect node_itr & node team_contact_ID
                   work_contacts_same_workplace_CS[team_contact_ID][1] = node_itr

                   # Connect team_contact_ID to other team members
                   team_contact_idx = 2 # Index for team member number connection corresponds to
                   for jj = 1:n_team_contacts
                       if ii != jj
                           jj_contact_ID = CS_workplace_contact_IDs[jj]
                           work_contacts_same_workplace_CS[team_contact_ID][team_contact_idx] = jj_contact_ID

                           # Increment index counter
                           team_contact_idx +=1
                       end
                   end

                   # Update assignment status for team member node
                   node_CS_assigned[team_contact_ID] = 1
                end

                # Update assignment status for node_itr
                node_CS_assigned[node_itr] = 1

            end
       end
   end

   # Create household contacts & assign a household ID
   csum_household_size_distribution = cumsum(household_size_distribution) # Get cumulative sum, used in each repeat of while loop
   worker_id = 1     # Initialise variables to be incremented in while loop
   household_id = 1
   while worker_id <= cmax
       household_size = findfirst(x->x>rand(rng), csum_household_size_distribution)  # generate random household size

       # check we don't exceed remaining workers left
       household_size = min(household_size, cmax-worker_id+1)

       # create fully-connected household
       for ii = 0:(household_size-1)
           for jj = 0:(household_size-1)
               # Get household contacts
               if ii != jj
                   push!(household_contacts[worker_id+ii], worker_id+jj)
                   household_contacts_per_node[worker_id+ii] += 1
               end
           end

           # Assign household ID to node ii
           worker_nodes[worker_id+ii].household_ID = household_id
       end

       worker_id += household_size
       household_id += 1
   end

   # Assign number of households in total to output variable
   n_households = household_id

   # Define what is returned from the function
   if CS_active_flag == true
       contacts = contacts_struct(cmax=cmax,endtime=endtime,work_contacts_same_workplace=work_contacts_same_workplace,
       work_contacts_other_workplace=work_contacts_other_workplace,
       social_contacts_per_node = social_contacts_per_node,
       work_contacts_same_workplace_CS = work_contacts_same_workplace_CS)

       return contacts::contacts_struct,
           social_contacts::Array{Array{Int64,1},1},
           household_contacts::Array{Array{Int64,1},1},
           work_contacts_same_workplace_per_node::Array{Int64,1},
           work_contacts_other_workplace_per_node::Array{Int64,1},
           work_contacts_same_workplace_per_node_CS::Array{Int64,1},
           household_contacts_per_node::Array{Int64,1},
           n_households::Int64
   else
       contacts = contacts_struct(cmax=cmax,endtime=endtime,work_contacts_same_workplace=work_contacts_same_workplace,
       work_contacts_other_workplace=work_contacts_other_workplace,
       social_contacts_per_node = social_contacts_per_node)

       return contacts::contacts_struct,
           social_contacts::Array{Array{Int64,1},1},
           household_contacts::Array{Array{Int64,1},1},
           work_contacts_same_workplace_per_node::Array{Int64,1},
           work_contacts_other_workplace_per_node::Array{Int64,1},
           household_contacts_per_node::Array{Int64,1},
           n_households::Int64
   end
end

# Premake dynamic worker contacts,
# to be loaded in ahead of simulation
function generate_dynamic_worker_contacts(RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            worker_nodes::Array{worker_params,1},
                                            dynamic_conts_mean::Array{Float64,1},
                                            dynamic_conts_sd::Array{Float64,1})
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# worker nodes - array with entry per worker
# dynamic_conts_mean, dynamic_conts_sd - Distribution properties for dynamic worker contacts

# Outputs:
# dynamic_worker_contacts - Per node, a record of dynamic worker contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_worker_contacts = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for node_itr = 1:cmax
        if worker_nodes[node_itr].returned_to_work==1 # Add dynamic links, if returned to work
            # Get dynamic worker group type for node_itr
            node_dynamic_grp_ID = worker_nodes[node_itr].sector_ID

            for time_itr = 1:endtime
                # Generate number of dynamic contacts from appropriate distribution
                gg = round(Int64,abs(randn(rng)*dynamic_conts_sd[node_dynamic_grp_ID]+dynamic_conts_mean[node_dynamic_grp_ID]))

                # If dynamic worker contacts made, assign to output variable
                if gg > 0
                    dynamic_worker_contacts[time_itr,node_itr] = zeros(gg)

                    # Generate required number of contacts
                    for contact_itr = 1:gg
                        gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        while gg1 == node_itr # Redraw if returned the index node
                            gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        end
                        dynamic_worker_contacts[time_itr,node_itr][contact_itr] = gg1
                    end
                else # No dynamic worker contacts made on given timestep
                    dynamic_worker_contacts[time_itr,node_itr] = Int64[]
                end
            end
        end
    end

    return dynamic_worker_contacts::Array{Array{Int64,1},2}

end

# Premake dynamic worker contacts,
# to be loaded in ahead of simulation
function generate_social_contacts_each_day(rng::MersenneTwister,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_contacts_per_node::Array{Int64,1},
                                            n_social_mean_workday::Int64,
                                            n_social_mean_nonworkday::Int64)
# Inputs:
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# social_contacts_per_node - Total amount of social contacts each individual has (entry per individual)
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those potentially have social contacts,
    assign a sample of those contacts for each timestep
    """
    for node_itr = 1:cmax

        # Add in social contacts if possible
        if social_contacts_per_node[node_itr] > 0

            # Social contacts may differ each day
            for time_itr = 1:endtime

                # Draw random number of social links to be made today
                n_social_workday = min(rand(rng,Poisson(n_social_mean_workday)), social_contacts_per_node[node_itr])
                n_social_nonworkday = min(rand(rng,Poisson(n_social_mean_nonworkday)), social_contacts_per_node[node_itr])

                if n_social_workday >0
                    workday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_workday)

                    for wd_social_it = 1:n_social_workday
                        wd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        workday_social_contacts_by_day[time_itr,node_itr][wd_social_it] = social_contacts[node_itr][wd_contact_idx]
                    end
                else
                    workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end

                if n_social_nonworkday >0
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_nonworkday)

                    for nwd_social_it = 1:n_social_nonworkday
                        nwd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        nonworkday_social_contacts_by_day[time_itr,node_itr][nwd_social_it] = social_contacts[node_itr][nwd_contact_idx]
                    end
                else
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end
            end
        else
            for time_itr = 1:endtime
                workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
            end
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


function generate_random_contacts(RNGseed::Int64,
                                    cmax::Int64,
                                    endtime::Int64,
                                    prob_random_contact::Float64)

    rng = MersenneTwister(RNGseed)

    random_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    for time_itr = 1:endtime
        for node_id1 = 1:(cmax-1)
            if isassigned(random_contacts_by_day, time_itr, node_id1) == false
                random_contacts_by_day[time_itr,node_id1] = Int64[]
            end
            for node_id2 = node_id1:cmax
                if isassigned(random_contacts_by_day, time_itr, node_id2) == false
                    random_contacts_by_day[time_itr,node_id2] = Int64[]
                end
                if rand(rng) < prob_random_contact
                    push!(random_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(random_contacts_by_day[time_itr,node_id2], node_id1)
                end
            end
        end
    end

    return random_contacts_by_day::Array{Array{Int64,1},2}
end

## Functions to regenerate specific layers of the network ##

# Were worker proportion to change from initial value to a non-zero value
# this function allows regeneration of returned_to_work for each node &
# reconstruct work networks.
function regenerate_worker_contacts(cmax::Int64,network_parameters::network_params,RNGseed::Int64)
# Inputs:
# cmax - Number of workers in the system
# network_parameters
#   worker nodes - array containing tuples of worker info
#   prob_workertype_contact - probability of making contact with others in worktype (DIFFERENT WORKPLACE)
#   prob_anyworker_contact - probability of making contact with others in DIFFERENT WORKTYPE
#   prob_social_contact - probability of making a contact socially
#   dd_within_workplace - mean degree of workers within workplace
#   household_size_distribution - distribution of household sizes
# RNGseed:Int64 - Value to seed the random number generator

# Outputs:
# work_contacts - Vector of vectors with IDs of contacts
# work_contacts_per_node - Total number of regular contacts within workplace

@unpack worker_nodes, workplace_sizes,prob_workertype_contact,prob_anyworker_contact,
    prob_social_contact,dd_within_workplace, household_size_distribution,
    workplace_info, CS_active_flag = network_parameters

    # Set the RNG
    rng = MersenneTwister(RNGseed)

    # Redo the return to workplace values
    for ii = 1:cmax
        returned_to_work = Int64(rand(rng)<workpercent[workertype_ids[ii]])   # decide if returning to work
        worker_nodes[ii].returned_to_work = returned_to_work  # add worker info to node
    end

    # Initialise vector of vectors storing IDs of contacts for each node
    # at workplace (non-covid secure setting)
    work_contacts_same_workplace = Array{Array{Int64,1},1}(undef,cmax)
    work_contacts_other_workplace = Array{Array{Int64,1},1}(undef,cmax)

    # If CS settings are active,
    # initialise vector of vectors storing IDs of contacts for each node
    # at workplace (covid-secure setting) & vector giving total contacts made by node
    if CS_active_flag == true
        work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
        work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
    end

    for ii = 1:cmax
        work_contacts_same_workplace[ii] = Int64[]
        work_contacts_other_workplace[ii] = Int64[]

        if CS_active_flag == true
            work_contacts_same_workplace_CS[ii] = Int64[]
        end
    end

    # Initialise vectors giving total contacts made by node
    work_contacts_same_workplace_per_node = zeros(Int64,cmax)
    work_contacts_other_workplace_per_node = zeros(Int64,cmax)

    # Construct network links for social, household, & workplace contacts
    # Workplace contacts done for non-covid secure setting
    for ii = 1:(cmax-1)

        returned_to_work::Int64 = worker_nodes[ii].returned_to_work
        worker_grp_idx::Int64 = worker_nodes[ii].sector_ID
        workplace_idx::Int64 = worker_nodes[ii].workplace_ID

        for jj = (ii+1):cmax

            ### On days when at work, increased contacts with others at work  ###
            if (returned_to_work == 1) & (worker_nodes[jj].returned_to_work==1)  # both returned to work

                ## For workers in the same group and workplace, edges form according to ER graph
                if (worker_nodes[jj].sector_ID == worker_grp_idx) & (worker_nodes[jj].workplace_ID == workplace_idx)
                    if rand(rng) < dd_within_workplace[worker_grp_idx]/(workplace_sizes[worker_grp_idx][workplace_idx] - 1) ## ER component
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_same_workplace[ii],jj)
                        push!(work_contacts_same_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_same_workplace_per_node[ii] += 1
                        work_contacts_same_workplace_per_node[jj] += 1
                    end

                ## same worker group, different workplace
            elseif (worker_nodes[jj].sector_ID == worker_grp_idx)

                    if rand(rng) < prob_workertype_contact[worker_grp_idx]
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_other_workplace[ii],jj)
                        push!(work_contacts_other_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_other_workplace_per_node[ii] += 1
                        work_contacts_other_workplace_per_node[jj] += 1
                    end

                ## different worker group
                else
                    if rand(rng) < prob_anyworker_contact
                        # Assign IDs of contact to tuple vectors
                        push!(work_contacts_other_workplace[ii],jj)
                        push!(work_contacts_other_workplace[jj],ii)

                        # Increment number of contacts for nodes ii & jj
                        work_contacts_other_workplace_per_node[ii] += 1
                        work_contacts_other_workplace_per_node[jj] += 1
                    end
                end
            end
        end
    end

    # If covid-secure workplaces in use,
    # construct in-workplace contacts assuming a covid secure setting
    if CS_active_flag == true
        node_CS_assigned = zeros(cmax)
        for node_itr = 1:(cmax-1) # Note, final node will have been accounted for

            # Check if node has already been assigned a workplace "team"
            # If not, go through assignment loop
             if node_CS_assigned[node_itr] == 0

                 # Get permitted team size
                 relevant_sector = worker_nodes[node_itr].sector_ID
                 relevant_workplace = worker_nodes[node_itr].workplace_ID
                 permitted_team_size = workplace_info[relevant_sector][relevant_workplace].CS_team_size

                 # Load pre-existing workplace contacts
                 nonCS_workplace_contacts = work_contacts_same_workplace[node_itr]

                 # Check how many workplaces have yet to be allocated
                 node_to_be_assigned = convert(Int64,sum(node_CS_assigned[nonCS_workplace_contacts].==0))
                 if node_to_be_assigned > 0
                     # Assign workplaces that could be allocated to vector
                     to_be_assigned_workplace_contacts = zeros(Int64,node_to_be_assigned)
                     assignment_idx = 1
                     for assign_itr = 1:length(nonCS_workplace_contacts)
                         ID_to_check = nonCS_workplace_contacts[assign_itr]
                         if node_CS_assigned[ID_to_check] == 0
                             to_be_assigned_workplace_contacts[assignment_idx] = ID_to_check

                             # Increment counter
                             assignment_idx += 1
                         end
                     end
                 else # No other workers available to be allocated to this team
                     to_be_assigned_workplace_contacts = Int64[]
                 end


                 # If number of existing contacts is below (max team size-1),
                 # create fully-connected team
                 n_to_be_assigned_workplace_contacts = length(to_be_assigned_workplace_contacts)

                 # If existing workplace contact amount larger than permitted size,
                 # pick a subset, and populate those vectors
                 if n_to_be_assigned_workplace_contacts > (permitted_team_size-1)
                     CS_workplace_contact_IDs = sample(rng, to_be_assigned_workplace_contacts, permitted_team_size-1; replace=false, ordered=false)
                 else
                     CS_workplace_contact_IDs = copy(to_be_assigned_workplace_contacts)
                 end
                 n_team_contacts = length(CS_workplace_contact_IDs)

                 # Populate workplace_CS contacts for node_itr
                 work_contacts_same_workplace_CS[node_itr] = CS_workplace_contact_IDs
                 work_contacts_same_workplace_per_node_CS[node_itr] = n_team_contacts


                 # Connect each team member to all other team members
                 for ii = 1:n_team_contacts

                    # Get ID of team member
                    team_contact_ID = CS_workplace_contact_IDs[ii]

                    # Initialise storage of workplace_CS contacts for node team_contact_ID
                    work_contacts_same_workplace_CS[team_contact_ID] = zeros(Int64,n_team_contacts)

                    # Set number of contacts for node_itr & node team_contact_ID
                    work_contacts_same_workplace_per_node_CS[team_contact_ID] = n_team_contacts

                    # Connect node_itr & node team_contact_ID
                    work_contacts_same_workplace_CS[team_contact_ID][1] = node_itr

                    # Connect team_contact_ID to other team members
                    team_contact_idx = 2 # Index for team member number connection corresponds to
                    for jj = 1:n_team_contacts
                        if ii != jj
                            jj_contact_ID = CS_workplace_contact_IDs[jj]
                            work_contacts_same_workplace_CS[team_contact_ID][team_contact_idx] = jj_contact_ID

                            # Increment index counter
                            team_contact_idx +=1
                        end
                    end

                    # Update assignment status for team member node
                    node_CS_assigned[team_contact_ID] = 1
                 end

                 # Update assignment status for node_itr
                 node_CS_assigned[node_itr] = 1

             end
        end
    end

    # Define what is returned from the function
    if CS_active_flag == true
        return  work_contacts_same_workplace::Array{Array{Int64,1},1},
                work_contacts_other_workplace::Array{Array{Int64,1},1},
                work_contacts_same_workplace_per_node::Array{Int64,1},
                work_contacts_other_workplace_per_node::Array{Int64,1},
                work_contacts_same_workplace_CS::Array{Array{Int64,1},1},
                work_contacts_same_workplace_per_node_CS::Array{Int64,1}
    else
        return  work_contacts_same_workplace::Array{Array{Int64,1},1},
                work_contacts_other_workplace::Array{Array{Int64,1},1},
                work_contacts_same_workplace_per_node::Array{Int64,1},
                work_contacts_other_workplace_per_node::Array{Int64,1}
    end
end

function regenerate_possible_social_contacts(workday_social_contacts_by_day::Array{Array{Int64,1},2},
                                            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2},
                                            prob_social_contact::Float64,
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            n_social_mean_workday::Int64,
                                            n_social_mean_nonworkday::Int64)
# Inputs:
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - Per node, a record of social contacts made on each day
# prob_social_contact - Connection probability with any other individual socially
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# starttime - Timestep to begin regenerating contacts from
# endtime - Number of timesteps simulation will run for
# social_contacts - array of arrays containing all possible social contacts for each individual
# n_social_mean_workday, n_social_mean_nonworkday - Distribution properties (mean value of Poisson) for social worker contacts

# Outputs:
# social_contacts_per_node - Total amount of social contacts each individual has (entry per individual)
# workday_social_contacts_by_day, nonworkday_social_contacts_by_day
#       - A now amended record of social contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Draw new collection of possible social contacts
    """
    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    social_contacts = Array{Array{Int64,1},1}(undef,cmax)

    # Initialise vectors giving total contacts made by node
    social_contacts_per_node = zeros(Int64,cmax)

    # Construct network links for social, household, & workplace contacts
    # Workplace contacts done for non-covid secure setting
    for ii = 1:(cmax-1)
        for jj = (ii+1):cmax

            ### Social contacts with anyone ###
            if rand(rng) < prob_social_contact
                # Assign IDs of contact to tuple vectors
                push!(social_contacts[ii],jj)
                push!(social_contacts[jj],ii)
                # Increment number of contacts for nodes ii & jj
                social_contacts_per_node[ii] += 1
                social_contacts_per_node[jj] += 1
            end
        end
    end

    """
    Iterate over all nodes
    For those potentially have social contacts,
    assign a sample of those contacts for each timestep
    """
    for node_itr = 1:cmax

        # Add in social contacts if possible
        if social_contacts_per_node[node_itr] > 0

            # Social contacts may differ each day
            for time_itr = starttime:endtime

                # Draw random number of social links to be made today
                n_social_workday = min(rand(rng,Poisson(n_social_mean_workday)), social_contacts_per_node[node_itr])
                n_social_nonworkday = min(rand(rng,Poisson(n_social_mean_nonworkday)), social_contacts_per_node[node_itr])

                if n_social_workday >0
                    workday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_workday)

                    for wd_social_it = 1:n_social_workday
                        wd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        workday_social_contacts_by_day[time_itr,node_itr][wd_social_it] = social_contacts[node_itr][wd_contact_idx]
                    end
                else
                    workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end

                if n_social_nonworkday >0
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_nonworkday)

                    for nwd_social_it = 1:n_social_nonworkday
                        nwd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        nonworkday_social_contacts_by_day[time_itr,node_itr][nwd_social_it] = social_contacts[node_itr][nwd_contact_idx]
                    end
                else
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end
            end

        end
    end

    return social_contacts_per_node::Array{Int64,1},
            workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


# Redraw dynamic worker contacts from specified start time to end time
# Allocate to pre-existing array
function regenerate_dynamic_worker_contacts(dynamic_worker_contacts::Array{Array{Int64,1},2},
                                            RNGseed::Int64,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            worker_nodes::Array{worker_params,1},
                                            dynamic_conts_mean::Array{Float64,1},
                                            dynamic_conts_sd::Array{Float64,1})

# Inputs:
# dynamic_worker_contacts - Per node, a record of dynamic worker contacts made on each day
# RNGseed - Seed the random number generator
# cmax - Number of nodes in the system
# starttime - Timestep to begin regenerating contacts from
# endtime - Number of timesteps simulation will run for
# worker nodes - array with entry per worker
# dynamic_conts_mean, dynamic_conts_sd - Distribution properties for dynamic worker contacts

# Outputs:
# dynamic_worker_contacts - The now amended record of dynamic worker contacts made on each day

    """
    Set the RNG
    """
    rng = MersenneTwister(RNGseed)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for node_itr = 1:cmax
        if worker_nodes[node_itr].returned_to_work==1 # Add dynamic links, if returned to work
            # Get dynamic worker group type for node_itr
            node_dynamic_grp_ID = worker_nodes[node_itr].sector_ID

            for time_itr = starttime:endtime
                # Generate number of dynamic contacts from appropriate distribution
                gg = round(Int64,abs(randn(rng)*dynamic_conts_sd[node_dynamic_grp_ID]+dynamic_conts_mean[node_dynamic_grp_ID]))

                # If dynamic worker contacts made, assign to output variable
                if gg > 0
                    dynamic_worker_contacts[time_itr,node_itr] = zeros(gg)

                    # Generate required number of contacts
                    for contact_itr = 1:gg
                        gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        while gg1 == node_itr # Redraw if returned the index node
                            gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                        end
                        dynamic_worker_contacts[time_itr,node_itr][contact_itr] = gg1
                    end
                else # No dynamic worker contacts made on given timestep
                    dynamic_worker_contacts[time_itr,node_itr] = Int64[]
                end
            end
        end
    end

    return dynamic_worker_contacts::Array{Array{Int64,1},2}
end