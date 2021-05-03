@everywhere using DrWatson
@everywhere @quickactivate

@everywhere include(srcdir("multiple_species.jl"))
#using Plots
@everywhere using AlertPushover, ProgressMeter



params = dict_list( Dict(
:h => [0.05, 0.1, 0.15, 0.2],
:ρ => [0.05, 0.1, 0.2, 0.3],
:T => [1.]
)) .|> (p->Multiparams(;pairs(p)...))


for param ∈ params
        #set initial dist
        initial_dist = centre_red(param)

        #set jump rates
        props_red       = Weights([0.,1., 1., 1., 1.])
        props_blue      = Weights([0.,1., 1., 1., 1.])
        rate_dictionary = Dict([("red",props_red),("blue",props_blue)])


        simparam = Simparams(param,initial_dist,rate_dictionary,false)

        sendto(workers(), simparam=simparam)
        @distributed for i in 1:100
                run_stepsave_const(simparam)
        end

        alert()
end


@alert @showprogress for param ∈ params
        initial_dist = centre_red(param)
        PDE_nodrift(param,initial_dist)
end



##
using JLD2
param = params[6]
@unpack h, ρ, Nr, Nb, T, τ, box = param
name = @ntuple  h ρ Nr Nb T τ box


y = fill(0.,param.width)
for i =1:96
        number = string(i)
        @load datadir("sims","multi_models_test", string(savename(name),"_#",number,".jld2")) model
        for i in 1:nagents(model)
                y[model[i].pos[1]] += 1
        end
end

@unpack Nr, Nb, h, Ω, timesteps, Δt, Dr, Db, T, saveon, box, dim, dr, db, α, ρ= param
name  = @ntuple h ρ Nr Nb T Δt box
@load datadir("PDE", savename(name, "jld2")) Pb


y = y/(param.Nr*96)


model.agents


##plot PDE
using Plots
plot()


y_data = Pr[1]
        y = fill(0.,param.width)
        for i = 1:param.width
                y[i]   = mean(y_data[i,:])
        end

plot!(range(0,stop = 10,length=param.width), y)

## plot model
using InteractiveDynamics
        import CairoMakie # choosing a plotting backend

        groupcolor(a) = a.colour == "red" ?  :red : :blue
        groupmarker(a) = a.colour == "red" ? :circle : :circle


        figure, _ = abm_plot(model; ac = groupcolor, am = groupmarker, as = 2, equalaspect= false)
        figure
