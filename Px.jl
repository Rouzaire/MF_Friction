@everywhere using DelimitedFiles , Distributed , Plotly , PyPlot, LaTeXStrings, StatsPlots
println() ; println() ; println()

@everywhere function PX_A(M)
    # Physical Parameters
        N   = 100 # number of particles
        mass   = 1 # mass of each particle (all blocks are identical for now)
        k_MF   = 1 # MF strength
        strength_potential = 1.0 # note that k_spring = 2 in our model
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 3*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    ## Creation of the distribution
        distribution_name = "gamma"
        k_gamma = 1.0
        distrib = [createCusps(distribution_name,k_gamma) for n in 1:N]

        filepath = ".\\PXA100\\"
        f_vel = open(filepath*"MeanVelocityContinuous_k"*string(k_gamma)*"_"*string(M)*".txt", "w")
        f_t = open(filepath*"effective_final_time_k"*string(k_gamma)*"_"*string(M)*".txt", "w")
        f_dtwell = open(filepath*"dtwell_k"*string(k_gamma)*"_"*string(M)*".txt", "w")

    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 1.5 # arbitrary, but it is bigger than the larger critical force encountered till now (11.5 for k = 1)
        t_transients = 50
        force = f_initial
        while t < t_transients
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_dtwell,string(getWellNumber(particle_position[1],distrib[1])) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : decrease force linearly at each timestep (approximation of continuity) until force = 0
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        df = f_initial*dt/(t_final - t_transients)
        under_threshold = fill(false,10)
        threshold = 0.05
        while t < t_final
            force = force - df
            t = t + dt
            write(f_dtwell,string(getWellNumber(particle_position[1],distrib[1])) * "\n")
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
            # Check every timestep for stopping condition
            # As soon as the mean velocity reaches 0 (in fact reaches the threshold N times), stop the simulation, print out effective final time and PX
            under_threshold = circshift(under_threshold, -1)
            under_threshold[end] = (getMeanVelocity(particle_velocity) < threshold)
            if under_threshold == fill(true,10)
                # print("Le break a fonctionné : ")
                break
            end # if
        end
        # println("on est sorti de la simulation.")
        # println("Effective Final Time = ",t)
        # println("Final Force = ",force)
        write(f_t,string(t) * "\n"*string(force) * "\n")

        close(f_vel)
        close(f_t)
        close(f_dtwell)

        PX_distribution = [getForceYield(particle_position[i],distrib[i]) for i in 1:N]
        PDelta_distribution = [getDelta(particle_position[i],distrib[i]) for i in 1:N]

        open(filepath*"PX_k"*string(k_gamma)*"_"*string(M)*".txt", "w") do io
               writedlm(io, [particle_position particle_velocity PX_distribution PDelta_distribution])
           end
end # function

@everywhere function PX_B(k_gamma,M)
    # Physical Parameters
        N   = 100 # number of particles
        mass   = 1 # mass of each particle (all blocks are identical for now)
        k_MF   = 1 # MF strength
        strength_potential = 1.0 # note that k_spring = 2 in our model
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 3*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    ## Creation of the distribution
        distribution_name = "gamma"
        distribution_param = 10/k_gamma # (param,k) = (10/x,x)
        distrib = [createCusps(distribution_name,distribution_param) for n in 1:N]

        filepath = ".\\PXB100\\"
        f_vel = open(filepath*"MeanVelocityContinuous_k"*string(k_gamma)*"_"*string(M)*".txt", "w")
        f_t = open(filepath*"effective_final_time_k"*string(k_gamma)*"_"*string(M)*".txt", "w")
        f_dtwell = open(filepath*"dtwell_k"*string(k_gamma)*"_"*string(M)*".txt", "w")

    # Phase 1 : wait for transients to fade out at high force
        force = 20
        t_transients = 50
        while t < t_transients
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_dtwell,string(getWellNumber(particle_position[1],distrib[1])) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : wait for transients to fade out at f = 0
        force = 0
        while t < 2*t_transients
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_dtwell,string(getWellNumber(particle_position[1],distrib[1])) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 3 : increase force linearly at each timestep (approximation of continuity) from force = 0 until stop
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        f_final = 20 # arbitrary, but it is bigger than the larger critical force encountered till now (11.5 for k = 1)
        df = f_final*dt/(t_final - 2*t_transients)
        above_threshold = fill(false,10)
        threshold = 2
        while t < t_final
            force = force + df
            t = t + dt
            write(f_dtwell,string(getWellNumber(particle_position[1],distrib[1])) * "\n")
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
            # Check every timestep for stopping condition
            # As soon as the mean velocity reaches 0 (in fact reaches the threshold N times), stop the simulation, print out effective final time and PX
            above_threshold = circshift(above_threshold, -1)
            above_threshold[end] = (getMeanVelocity(particle_velocity) > threshold)
            if above_threshold == fill(true,10)
                # print("Le break a fonctionné : ")
                break
            end # if
        end
        # println("on est sorti de la simulation.")
        print("Effective Final Time = ",t)
        println("     Final Force = ",force)
        write(f_t,string(t) * "\n"*string(force) * "\n")

        close(f_vel)
        close(f_t)
        close(f_dtwell)

        PX_distribution = [getForceYield(particle_position[i],distrib[i]) for i in 1:N]
        PDelta_distribution = [getDelta(particle_position[i],distrib[i]) for i in 1:N]

        open(filepath*"PXB_k"*string(k_gamma)*"_"*string(M)*".txt", "w") do io
               writedlm(io, [particle_position particle_velocity PX_distribution PDelta_distribution])
           end
end # function

filepath = ".\\PXA100\\"
k_gamma = 1.0
M = 100
@time pmap(PX_A,[i for i in 101:200])

px_A = pdelta_A = []
for m in 1:M
    global px_A = vcat(px_A,readdlm(filepath*"PX_k"*string(k_gamma)*"_"*string(m)*".txt", Float64)[:,3])
    global pdelta_A = vcat(pdelta_A,readdlm(filepath*"PX_k"*string(k_gamma)*"_"*string(m)*".txt", Float64)[:,4])
end
pyplot()
Plots.histogram(px_A,xlabel = "\$x_{\\sigma}\$",ylabel = "\$P (x_{\\sigma})\$")
Plots.plot(StatsPlots.fit(Gamma,px_A))

# Plots.histogram(pdelta_A,normalize=:pdf ,nbin = 100,label="Width of the stopping wells")
# Plots.histogram!(0.5*createCusps("gamma",5),normalize=:pdf ,nbin=50,label="Gamma Distribution")
# xlabel!("\$\\Delta\$")
# ylabel!("\$P(\\Delta)\$")
# Plots.pdf("Stopping Delta Distribution")


wellnb = readdlm(filepath*"dtwell_k"*string(k_gamma)*"_"*string(1)*".txt")
d_well = [wellnb[i+1] - wellnb[i] for i in 1:length(wellnb)-1]
Plots.histogram(d_well,legend=nothing,yaxis = (:log))
Plots.pdf("Distribution_D_well")

final_force_k1 = fill(-1.0,M)
# final_force_k2 = fill(-1.0,M)
for m in 1:M
    final_force_k1[m] = readdlm(".\\PXA100\\effective_final_time_k"*string(1.0)*"_"*string(m)*".txt", Float64)[2]
    # final_force_k2[m] = readdlm(".\\PX100\\effective_final_time_k"*string(2.0)*"_"*string(m)*".txt", Float64)[2]
end
Plots.histogram(final_force_k1,nbins=25,normalize=:pdf,label="k=1, N = 100")
# Plots.histogram!(final_force_k2,normalize=:pdf,nbins=15,label="k=2, N = 100")
title!("Distribution of force when simulation stops")
xlabel!(L"f")
ylabel!(L"P(f)")
Plots.pdf("Distribution_force_stop_1_2")
