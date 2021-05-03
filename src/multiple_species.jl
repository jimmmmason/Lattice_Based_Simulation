using DrWatson
@quickactivate "Exclusion Process"

using Agents, LinearAlgebra, Random, StatsBase, Parameters, ProgressMeter


## Parameters

@with_kw struct Multiparams @deftype Float64
    box::Tuple{Float64,Float64} = (10,2)
    h  = 0.1
    ρ  = 0.2
    T  = 1
    Δt = 0.01
    Dr = 1
    Db = 1
    saveon::Bool = true

    width ::Int64 = round(box[1]/h)
    height::Int64 = round(box[2]/h)
    dim ::Tuple{Int64,Int64}= (width, height)

    N = 20*ρ/h^2
    Nr::Int64 = round(0.2*N)
    Nb::Int64 = round(0.8*N)
    Ω::Array{Array{Int64,1},2} = [[i,j] for j ∈ 1:height, i ∈ 1:width]
    α = π/2

    timesteps::Int64 =  round(T/Δt)

    τ      = h^2/(Nr*Dr+Nb*Db)
    chunks::Int64 = 10
    nhops ::Int64  = round(T/(τ*chunks))
    dr = 2*Dr/(Dr+Db)
    db = 2*Db/(Dr+Db)
end

@with_kw struct Simparams @deftype Float64
        params::Multiparams
        initial_dist::Dict
        rate_dictionary::Dict
        overlap::Bool = false
end



## Simulation

#set up
mutable struct MultipleSpeciesAgent <: AbstractAgent
        id::Int # The identifier number of the agent
        pos::Dims{2} # The x, y location of the agent on a 2D grid
        colour::String #total displacement over time
end


#Model step picks the particels in a random order and applies agent_step

#set up for a constant rate function
function model_step_const!(model)
        agent        = random_agent(model)
        jump_weights = Weights(model.rates[agent.colour])
        jump         = sample([(0,0),(1,0), (0,1), (-1,0), (0,-1)],jump_weights)
        walk!(agent, jump, model; ifempty = true)
        model.tick  += 1
end

#set up for a dynamic rate function
function model_step_dyn!(model)
        agent        = random_agent(model)
        jump_weights = Weights(model.rates[agent.colour](collect(agent.pos)))
        jump         = sample([(0,0),(1,0), (0,1), (-1,0), (0,-1)],jump_weights)
        walk!(agent, jump, model; ifempty = true)
        model.tick  += 1
end
#initial dist with box of red in the middle
function centre_red(params)
        @unpack h, dim = params

        lower_mid = Int(round(4/h))+1
        upper_mid = Int(round(6/h))
        middle    = lower_mid:1:upper_mid
        dim_mid   = (Int(round(2/h)),Int(round(2/h)))

        Pr = fill(0.,dim)
        Pr[middle ,:] = fill(1/4,dim_mid)
        Pb = fill(1/16,dim)
        Pb[middle ,:] = fill(0.,dim_mid)

        return Dict(:"red"=> Pr, :"blue"=> Pb)
end

#initialises model
function multi_initialise(params::Multiparams, initial_dist::Dict, rate_dictionary::Dict, overlap = true)
    #dist is an array contianing 2 weight vectors
    @unpack height, width, dim, Nr, Nb, Ω = params

    #environment
    properties = @dict
    properties[:tick]  = 0
    properties[:rates] = rate_dictionary

    space = GridSpace(dim, periodic = true)
    model = ABM(MultipleSpeciesAgent, space; properties=properties)

    # Initial state
    Ω = [(i,j) for j in 1:height for i in 1:width]


    initial_distribution_r = sample(Ω, Weights(vec(initial_dist["red"])), Nr; replace=false, ordered=true)


    # if Pb doesn't overlap Pr we can reset it so not to missorder Pb entries
    if overlap == false
        Ω = [(i,j) for j in 1:height for i in 1:width]
    end

    initial_distribution_b = sample(Ω, Weights(vec(initial_dist["blue"])), Nb; replace=false, ordered=true)



    # add agents to model
    id   = 1 ::Int
    for x in initial_distribution_r
            agent = MultipleSpeciesAgent(id, x, "red")
            add_agent_pos!(agent, model)
            id+=1
    end

    for x in initial_distribution_b
            agent = MultipleSpeciesAgent(id, x, "blue")
            add_agent_pos!(agent, model)
            id+=1
    end

    return model
end

#runs saving every nhops for chunks
function run_stepsave_const(simparams::Simparams)
        @unpack params, initial_dist, rate_dictionary, overlap = simparam
        model = multi_initialise(params, initial_dist, rate_dictionary, overlap)
        @unpack nhops, box, chunks, h, ρ, Nr, Nb,  τ, box, saveon =params
        for i ∈ 1:chunks
                run!(model,dummystep,model_step_const!,nhops)
                if saveon
                        local T, ticks
                        ticks = model.tick
                        T = τ*ticks
                        name = @ntuple  h ρ Nr Nb T τ box
                        data       = Dict([("model",model)])
                        safesave(datadir("sims","multi_models_test", savename(name, "jld2")),data)
                end
        end
        return model
end



## Difference eq. solver

function neighbour(x,dim,Ω)
        neighbours = []
        for y ∈ Ω
            if norm(x-y)==1
                push!(neighbours,y)
            elseif (abs(x[1]-y[1]) == dim[1]-1)&(x[2]==y[2]) #i.e x is on the boundary
                push!(neighbours,y)
            elseif (abs(x[2]-y[2]) == dim[2]-1)&(x[1]==y[1])
                push!(neighbours,y)
            end
        end
        return neighbours
end


#need to chnage f into a vector I was seeing it as 2d for some reason...
function PDE_nodrift(params::Multiparams,initial_dist::Dict)
            @unpack Nr, Nb, h, Ω, timesteps, Δt, Dr, Db, T, saveon, box, dim, dr, db, α, ρ= params
            Pr = initial_dist["red"]
            Pb = initial_dist["blue"]
            @showprogress for step ∈ 1:timesteps
                local r,b
                r::Array{Float64,2}= fill(0.,dim)
                b::Array{Float64,2} = fill(0.,dim)
                for x ∈ Ω, y ∈ neighbour(x,dim,Ω)
                        local laplacer, EPr, laplaceb, EPb, Er, Eb
                        #one particle operator
                        laplacer::Float64 = Dr*(1/h^2)*(Pr[y...]-Pr[x...]) #λr(x,y)*(Pr[y...]-Pr[x...])  in the general case
                        #Interaction term
                        #gradient terms in E (leave h for now)
                        EPr::Float64      = (-1-dr*α)*(1/2)*(Pb[y...]+Pb[x...])*(Pr[y...]-Pr[x...])  +  (1+db*α)*(1/2)*(Pr[y...]+Pr[x...])*(Pb[y...]-Pb[x...])
                        #opposing colour drift term
                        #Efr      = (1+α*dr)*sgn(y-x)*fr*(1/2)*(Pr[y...]*Pb[y...]+Pr[x...]*Pb[x...])  +  (-α*db)*sgn(y-x)*fb*(1/2)*(Pr[y...]*Pb[y...]+Pr[x...]*Pb[x...])
                        #mean field self drift term
                        #Err      = sgn(y-x)*fr*(1/2)*(Pr[y...]*Pr[y...]+Pr[x...]*Pr[x...])
                        #full interaction term, correcting for h
                        #Er       = Dr*(EPr+Nb*(EPr+h*Efr)+(Nr-1)*h*Err)
                        Er::Float64      = Dr*Nb*EPr


                        laplaceb::Float64 = (Db/h^2)*(Pb[y...]-Pb[x...]) #λb(x,y)*(Pb[y...]-Pb[x...])  in the genebal case
                        EPb::Float64      = (-1-db*α)*(1/2)*(Pr[y...]+Pr[x...])*(Pb[y...]-Pb[x...])+(1+dr*α)*(1/2)*(Pb[y...]+Pb[x...])*(Pr[y...]-Pr[x...])
                        #Efb      = (1+α*db)*sgn(y-x)*fb*(1/2)*(Pb[y...]*Pr[y...]+Pb[x...]*Pr[x...])+(-α*dr)*sgn(y-x)*fr*(1/2)*(Pb[y...]*Pr[y...]+Pb[x...]*Pr[x...])
                        #Ebb      = sgn(y-x)*fb*(1/2)*(Pb[y...]*Pb[y...]+Pb[x...]*Pb[x...])
                        #Eb       = Db*(EPb+Nr*(EPb+h*Efb)+(Nb-1)*h*Ebb)
                        Eb::Float64       = Dr*Nr*EPb

                        r[x...] += Δt*(laplacer+Er)
                        b[x...] += Δt*(laplaceb+Eb)
                end

                Pr += r
                Pb += b

            end

            if saveon
                    data  = Dict([("Pr",Pr),("Pb",Pb)])
                    name  = @ntuple h ρ Nr Nb T Δt box
                    safesave(datadir("PDE", savename(name, "jld2")),data)
            end
            return [Pr,Pb]
        end



## Random
function normish(x,sig)
    v = fill(0.,length(x))
    for i in 1:length(x)
        local w
        w = exp(-norm(collect(x[i])-[height, width]./2)^2/(2*sig))
        v[i] = w
    end
  return v
end

function sgn(x)
    for i ∈ 1:length(x)
        if x[i] != 0
            return x[i]
        end
    end
end
