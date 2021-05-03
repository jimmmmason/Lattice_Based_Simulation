using DrWatson
@quickactivate "Exclusion Process"
using Distributed

machines = ["squiffy","coot",
            "clouds","cardioid","krisskross","diogenes"]

machines = ["nitrogen","cerium","statuepark","wheatear",
            "pluto","wagner"]

machines =  ["epicurus","verdi",
            "chopin","swan","ovenbird","tangram"]

machines = [ "bluetail","squiffy","citril","crow",
               "cuckoo","dotterel","gargoyle","garnish"]

addprocs(["garnish"])
addprocs(machines)
addprocs([("oblong",:auto)])

@sync @distributed for i = 1:nworkers()
                     host = gethostname()
                     println(host)
                  end




Pkg.add("ParallelDataTransfer")
@everywhere using ParallelDataTransfer
















@everywhere script = "/store/DAMTP/jm2386/Lattice_Based_Simulation/scripts/test_script.jl"
@everywhere include(script)

@distributed for i in 1:10
                println(i)
             end

rmprocs(procs()[2]:procs()[length(machines)+1])
rmprocs(1:18)

##



for height ∈ heights, width ∈ widths#, rho ∈ rhos
       @distributed  for (batch,i) in collect(Iterators.product(1:batches, 1:5))
                                local rho
                                rho = 1-0.01*i
                                println(rho)
                                #model = tagged_initialise(height,width,rho)
                                #ticks   = Int(100*rho*height*width) #so that t is proportional to 100
                                #run!(model,dummystep,model_step!,ticks)
                                #diff = diffusion_tensor_estimate(model)
                                #simulation = @ntuple batch ticks rho height width
                                #data       = Dict([("Diffusion_estimate",diff)])

                                #safesave(datadir("sims", "tagged", savename(simulation, "jld2")),data)
                                println("I contributed")
                                rho =5
                                println(rho)
                                println(batch)
         end
end



@sync @distributed for i = 1:nworkers()
                     host = gethostname()
                     println(host)
                  end
