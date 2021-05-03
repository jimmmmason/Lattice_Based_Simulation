using DrWatson
@quickactivate "Exclusion Process"
include(srcdir("multiple_species.jl"))
using StatsBase

τ     = 0.001
        dt = 0.01
        hs    = [0.05, 0.1, 0.125, 0.15, 0.175, 0.2]
        rhos  = [0.05, 0.1, 0.2, 0.3]

        #constant
        α = π/2

        #Diffusion parameters

        Dr = 1
        Db = 1
        fr = 0
        fb = 0

@alert @distributed  for h in hs
        for ρ in rhos
                N = 20*ρ/h^2
                Nr = Int(round(0.2*N))
                Nb = Int(round(0.8*N))

                width  = Int(round(10/h))
                height = Int(round(2/h))
                dim   = (width,height) #i.e real 'space' is (10,2)
                Ω     = [[i,j] for i ∈ 1:dim[1], j ∈ 1:dim[2]]


                #initial condidtions

                lower_mid = Int(round(4/h))+1
                upper_mid = Int(round(6/h))
                midd    = lower_mid:1:upper_mid
                dim_mid   = (Int(round(2/h)),Int(round(2/h)))

                Pr = fill(0.,dim)
                Pr[midd ,:] = fill(1/4,dim_mid)
                Pb = fill(1/16,dim)
                Pb[midd ,:] = fill(0.,dim_mid)

                props_red  = dt.*fill([0.,1., 1., 1., 1.],dim)
                props_blue = dt.*fill([0.,1., 1., 1., 1.],dim)

                rates = Dict([("red",props_red),("blue",props_blue)])
                #in future define function rates that takes agent as input

                Ω = [(i,j)  for j in 1:height for i in 1:width]
                dist = [Weights(reshape(Pr,length(Pr))), Weights(reshape(Pb,length(Pb)))]

                model = multi_initialise_2(height,width,Nr,Nb,dist,rates, false)

                for i in 1:100
                        run!(model,dummystep,model_step!,100)

                        ticks = model.tick
                        simulation = @ntuple ticks rho height width
                        data       = Dict([("Diffusion_estimate",diff)])

                        safesave(datadir("sims", "multi_box_10x2", savename(simulation, "jld2")),data)
                end
        end
end



using InteractiveDynamics
import CairoMakie # choosing a plotting backend

groupcolor(a) = a.colour == "red" ?  :red : :blue
groupmarker(a) = a.colour == "red" ? :circle : :circle


figure, _ = abm_plot(model; ac = groupcolor, am = groupmarker, as = 2, equalaspect= false)
figure

abm_video(
    "/store/DAMTP/jm2386/Lattice_Based_Simulation/plots/red_blue.mp4", model, dummystep, model_step!;
    ac = groupcolor, am = groupmarker, as = 2,
    framerate = 100, frames = 100000,
    title = "Schelling's segregation model"
)

using Plots

figureaxisplot = scatter(rand(100, 2))
figure = figureaxisplot.figure
