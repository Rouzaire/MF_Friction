@everywhere using DelimitedFiles , Statistics, Distributed , Plots, LaTeXStrings, StatsPlots, StatsBase, JLD
println() ; println() ; println()

@everywhere function Hysteresis_Above_Single_Realisation(distrib, mass, N,k_gamma,k_MF)
    # Physical Parameters
        strength_potential = 1.0
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)

        filepath = ".\\hystA\\"
        f_vel = open(filepath*"MeanVelocity_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*string(k_MF)*".txt", "w")
        f_f = open(filepath*"force_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*string(k_MF)*".txt", "w")


    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 2 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : decrease force linearly at each timestep (approximation of continuity) until force = 0
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        df = f_initial*dt/(t_final - t_transients)
        while t < t_final
            force = force - df
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_f,string(force) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

        close(f_vel)
        close(f_f)
end # function
@everywhere function Hysteresis_Below_Single_Realisation(distrib, mass, N,k_gamma,k_MF)
    # Physical Parameters
        strength_potential = 1.0
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)

        filepath = ".\\hystB\\"
        f_vel = open(filepath*"MeanVelocity_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*string(k_MF)*".txt", "w")
        f_f = open(filepath*"force_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*string(k_MF)*".txt", "w")

    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 1.5 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : wait for transients to fade out at f = 0
        force = 0 # arbitrary, but it is bigger than the larger critical force encountered till now
        while t < 2*t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 3 : increase force linearly at each timestep (approximation of continuity) from force = 0 until stop
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        f_final = 2 # arbitrary, but it is bigger than the larger critical force (for below) encountered till now
        df = f_final*dt/(t_final - 2*t_transients)
        while t < t_final
            force = force + df
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_f,string(force) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

        close(f_vel)
        close(f_f)
end # function
@everywhere function Hysteresis_Above_Average_Disorder(N,simulation_Number)
    # Physical Parameters
        mass = 1
        k_MF = 1
        strength_potential = 1.0
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    # Creation of the distribution
        name = "gamma"
        k_gamma = 1
        distrib = [createCusps(name,k_gamma) for n in 1:N]


        filepath = ".\\hystAVG\\"
        f_vel = open(filepath*"MeanVelocity_A_"*simulation_Number*".txt", "w")
        f_f = open(filepath*"force_A_"*simulation_Number*".txt", "w")


    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 2 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : decrease force linearly at each timestep (approximation of continuity) until force = 0
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        df = f_initial*dt/(t_final - t_transients)
        while t < t_final
            force = force - df
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_f,string(force) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

        close(f_vel)
        close(f_f)
end # function
@everywhere function Hysteresis_Below_Average_Disorder(N,simulation_Number)
    # Physical Parameters
        mass = 1
        k_MF = 1
        strength_potential = 1.0
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    # Creation of the distribution
        name = "gamma"
        k_gamma = 1
        distrib = [createCusps(name,k_gamma) for n in 1:N]



        filepath = ".\\hystAVG\\"
        f_vel = open(filepath*"MeanVelocity_B_"*simulation_Number*".txt", "w")
        f_f = open(filepath*"force_B_"*simulation_Number*".txt", "w")



    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 1.5 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : wait for transients to fade out at f = 0
        force = 0 # arbitrary, but it is bigger than the larger critical force encountered till now
        while t < 2*t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 3 : increase force linearly at each timestep (approximation of continuity) from force = 0 until stop
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        f_final = 2 # arbitrary, but it is bigger than the larger critical force (for below) encountered till now
        df = f_final*dt/(t_final - 2*t_transients)
        while t < t_final
            force = force + df
            t = t + dt
            write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            write(f_f,string(force) * "\n")
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

        close(f_vel)
        close(f_f)
end # function
@everywhere function Measurements_Stop_Above(k_MF, k_gamma, ID)
    # Physical Parameters
        N = 100
        strength_potential = 1.0
        eta = 0.34641016151377546
        mass = 1000
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    ## Creation of the distribution
        distribution_name = "gamma"
        distrib = [createCusps(distribution_name,k_gamma) for n in 1:N]
        filepath = ".\\PXscan\\scan_mass\\"
        # f_vel = open(filepath*"MeanVelocity_k"*string(k_gamma)*"_m"*string(mass)*".txt", "w")
        # f_f = open(filepath*"force_k"*string(k_gamma)*"_m"*string(mass)*".txt", "w")

    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 2 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
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
            # write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
            # write(f_f,string(force) * "\n")
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
            end
        end
        PX_distribution = [getForceYield(particle_position[i],distrib[i]) for i in 1:N]
        PDelta_distribution = [getDelta(particle_position[i],distrib[i]) for i in 1:N]
        open(filepath*"PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N"*string(N)*"m"*string(mass)*"_"*ID*".txt", "w") do io
               writedlm(io, [particle_position particle_velocity PX_distribution PDelta_distribution])
        end # ends do block
        open(filepath*"stop_force"*string(k_gamma)*"_MF"*string(k_MF)*"_N"*string(N)*"m"*string(mass)*"_"*ID*".txt", "w") do io
               writedlm(io, [force t])
        end # ends do block
end # function
@everywhere function Measurements_Stop_Below(mass, N, ID)
    # Physical Parameters
        k_MF   = 1 # MF strength
        strength_potential = 1.0
        eta = 0.34641016151377546
    # Numerical Parameters
        dt = .1
        t_final = 2*1e4
        t = 0.0
        particle_position = fill(1e-4,N)
        particle_velocity = zeros(N)
    ## Creation of the distribution
        distribution_name = "gamma"
        k_gamma = 1.0
        distrib = [createCusps(distribution_name,k_gamma) for n in 1:N]

        filepath = ".\\StopB\\"

    # Phase 1 : wait for transients to fade out at f = f_initial
        f_initial = 2 # arbitrary, but it is bigger than the larger critical force encountered till now
        t_transients = 200
        force = f_initial
        while t < t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 2 : wait for transients to fade out at f = 0
        force = 0 # arbitrary, but it is bigger than the larger critical force encountered till now
        while t < 2*t_transients
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
        end

    # Phase 3 : increase force linearly at each timestep (approximation of continuity) until force = 2
        # df/dt should be sufficiently small so that it's quasi stationnary (and therefore transients effects should be negligeable in this Phase 2)
        f_final = 2 # arbitrary, but it is bigger than the larger critical force encountered till now
        df = f_final*dt/(t_final - 2*t_transients)
        above_threshold = fill(false,10)
        threshold = 2
        while t < t_final
            force = force + df
            t = t + dt
            for j in 1:N
                particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
            end
            # Check every timestep for stopping condition
            # As soon as the mean velocity departs from 0 (in fact reaches the threshold N times), stop the simulation, print out effective final time and PX
            above_threshold = circshift(above_threshold, -1)
            above_threshold[end] = (getMeanVelocity(particle_velocity) > threshold)
            if above_threshold == fill(true,10)
                break
            end
        end
        PX_distribution = [getForceYield(particle_position[i],distrib[i]) for i in 1:N]
        PDelta_distribution = [getDelta(particle_position[i],distrib[i]) for i in 1:N]
        open(filepath*"PX_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_"*ID*".txt", "w") do io
               writedlm(io, [particle_position particle_velocity PX_distribution PDelta_distribution])
        end # ends do block
        open(filepath*"stop_force"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_"*ID*".txt", "w") do io
               writedlm(io, [force t])
        end # ends do block
end # function

# N = 100
# M = 1000
# @time pmap(Hysteresis_Above_Average_Disorder,fill(N,M),[string(m) for m in 1:M])
# @time pmap(Hysteresis_Below_Average_Disorder,fill(N,M),[string(m) for m in 1:M])

# Creation of the distribution
# N = 100
# distribution_name = "gamma"
# k_gamma = 1.0 # k_gamma
# mass = 1
# distrib = [createCusps(distribution_name,k_gamma) for n in 1:N]

# k_MF = 10
# @time Hysteresis_Above(distrib,mass,N,k_gamma,k_MF)
# @time Hysteresis_Below(distrib,mass,N,k_gamma,k_MF)
#
# k_MF = 0.1
# @time Hysteresis_Above(distrib,mass,N,k_gamma,k_MF)
# @time Hysteresis_Below(distrib,mass,N,k_gamma,k_MF)
#
# k_MF = 1
# @time Hysteresis_Above(distrib,mass,N,k_gamma,k_MF)
# @time Hysteresis_Below(distrib,mass,N,k_gamma,k_MF)

# Distribution of Stopping Forces and P(x)
M = 1000
k_MF = 2
k_gamma = 1.0
@time pmap(Measurements_Stop_Above,fill(k_MF,M),fill(k_gamma,M),[string(i) for i in 1:M])
filepath = ".\\PXscan\\scan_mass\\"
PXA10 = Float64[]
PXA1000 = Float64[]
PXA1 = Float64[]
PXA01 = Float64[]
for i in 1:M
    global PXA10 = vcat(PXA10,readdlm(filepath*"PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N100m10_"*string(i)*".txt")[:,3])
    global PXA1000 = vcat(PXA1000,readdlm(filepath*"PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N100m1000_"*string(i)*".txt")[:,3])
    global PXA1 = vcat(PXA1,readdlm(filepath*"PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N100m1_"*string(i)*".txt")[:,3])
    global PXA01 = vcat(PXA01,readdlm(filepath*"PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N100m0.1_"*string(i)*".txt")[:,3])
end
# open(filepath*"m1N100.txt", "w") do io
#            writedlm(io, [PXA])
# end
Plots.histogram(PXA1,nbins=1000)
h = StatsBase.fit(Histogram,log10.(PXA1),nbins=1000)
x,y = [collect(h.edges[1]), h.weights]
Plots.plot(x[1:end-1],log10.(y),label="m = 1000b")
Plots.plot(x[1:end-1],log10.(y) ./ x[1:end-1],label="m = 100b",xlims=((-3,-1)))
ylims!((-2,0))
## Data analysis for Hysteresis curves (velocity vs force) for a single realisation
    # Plots.plot()
    # for kmeanfield in ["10", "1", "0.1"]
    #     println(kmeanfield)
    #     velA = readdlm(".\\hystA\\MeanVelocity_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*kmeanfield*".txt")
    #     velB = readdlm(".\\hystB\\MeanVelocity_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*kmeanfield*".txt")
    #     fA = readdlm(".\\hystA\\force_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*kmeanfield*".txt")
    #     fB = readdlm(".\\hystB\\force_k"*string(k_gamma)*"_m"*string(mass)*"_N"*string(N)*"_MF"*kmeanfield*".txt")
    #     display(Plots.plot!(fA,velA,label="Above, k_{MF} = "*kmeanfield))
    #     display(Plots.plot!(fB,velB,label="Below, k_{MF} = "*kmeanfield))
    # end
    # xlabel!("f")
    # ylabel!(L"\overline{v}(f)")
    # Plots.pdf("Hysteresis_diff_kMF_k"*string(Int(k_gamma))*"_N"*string(N)*".pdf")
    # Plots.plot()
    # for masse in ["10", "1", "0.1"]
    #     println(masse)
    #     velA = readdlm(".\\hystA\\MeanVelocity_k"*string(k_gamma)*"_m"*masse*"_N"*string(N)*".txt")
    #     velB = readdlm(".\\hystB\\MeanVelocity_k"*string(k_gamma)*"_m"*masse*"_N"*string(N)*".txt")
    #     fA = readdlm(".\\hystA\\force_k"*string(k_gamma)*"_m"*masse*"_N"*string(N)*".txt")
    #     fB = readdlm(".\\hystB\\force_k"*string(k_gamma)*"_m"*masse*"_N"*string(N)*".txt")
    #     display(Plots.plot!(fA,velA,label="Above, m = "*masse))
    #     display(Plots.plot!(fB,velB,label="Below, m = "*masse))
    # end
    # xlabel!("f")
    # ylabel!(L"\overline{v}(f)")
    # Plots.pdf("Hysteresis_k"*string(Int(k_gamma))*"_N"*string(N)*".pdf")

## Data analysis for Hysteresis curves (velocity vs force) averaged over disorder
    # forceA = readdlm(".\\hystAVG\\force_A_1.txt") ; forceB = readdlm(".\\hystAVG\\force_B_1.txt") # all force files are identical since the protocol is always the same
    # velA = zeros(length(forceA),M) ; velB = zeros(length(forceB),M)
    # for m in 1:M
    # velA[:,m] = readdlm(".\\hystAVG\\MeanVelocity_A_"*string(m)*".txt")
    # velB[:,m] = readdlm(".\\hystAVG\\MeanVelocity_B_"*string(m)*".txt")
    # end
    # meanA = mean(velA,dims=2) ; meanB = mean(velB,dims=2)
    # stdA = Statistics.std(velA,dims=2) ; stdB = Statistics.std(velB,dims=2)
    # minA = [abs.(minimum(velA[i,:])-meanA[i]) for i in 1:length(forceA)] ; maxA = [abs.(maximum(velA[i,:])-meanA[i]) for i in 1:length(forceA)]
    # minB = [abs.(minimum(velB[i,:])-meanB[i]) for i in 1:length(forceB)] ; maxB = [abs.(maximum(velB[i,:])-meanB[i]) for i in 1:length(forceB)]
    # # Computes ratio stopped
    # ratioA = zeros(length(forceA))
    # ratioB = zeros(length(forceB))
    # for i in 1:length(forceA)
    # tmp = velA[i,:]
    # ratioA[i] = length(tmp[tmp .< 0.01])/M*100
    # end
    # reverse!(ratioA)
    # for i in 1:length(forceB)
    # tmp = velB[i,:]
    # ratioB[i] = length(tmp[tmp .< 0.01])/M*100
    # end
    #
    #
    # pyplot()
    # increment = 3000
    # Plots.plot(box=true,legend=:left)
    # Plots.plot!(forceA[1:increment:end],meanA[1:increment:end],yerr=(minA[1:increment:end],maxA[1:increment:end]),label="Above",ribbon=stdA[1:increment:end],fillalpha=.5)
    # Plots.plot!(1:2, NaN.*(1:2), label = "Stop Ratio Above [%]", linecolor=:blue, ls=:dash, grid=false)
    # Plots.plot!(Plots.twinx(),forceA[1:Int(increment/10):end],ratioA[1:Int(increment/10):end],ls=:dash,lc=:blue,legend=false,ymirror=true)
    # xlabel!(L"f")
    # ylabel!(L"\langle \ \overline{\,v \,}\ \rangle_{disorder}")
    # # Plots.pdf("HystAVG_Above")
    #
    #
    # Plots.plot(box=true,legend=:left)
    # Plots.plot!(forceA[1:100:end],meanA[1:100:end],label="Above")
    # Plots.plot!(forceB[1:increment:end],meanB[1:increment:end],yerr=(minB[1:increment:end],maxB[1:increment:end]),label="Below",color=RGBA(0.8888735002725198,0.43564919034818983,0.2781229361419438,1.0),ribbon=stdB[1:increment:end],fillalpha=.5)
    # yaxis!((-0.16,6))
    # Plots.plot!(1:2, NaN.*(1:2), label = "Stop Ratio Below [%]", linecolor=:red, ls=:dash, grid=false)
    # Plots.plot!(Plots.twinx(),forceB[1:Int(increment/10):end],ratioB[1:Int(increment/10):end],ls=:dash,lc=:red,legend=false,ymirror=true)
    # xlabel!(L"f")
    # ylabel!(L"\langle \ \overline{\,v \,}\ \rangle_{disorder}")
    # Plots.pdf("HystAVG_Below")

# Justification of the choice of N : convN
MM = [500, 1000, 667, 500, 250, 166, 105]
NN = string.([20,50, 75, 100, 200, 300, 500])
Plots.plot()
for n in length(MM):-1:1
    PXA = Float64[]
    for i in 1:MM[n]
        PXA = vcat(PXA,readdlm(".\\convN\\PX_k"*string(1.0)*"_MF"*string(2)*"_N"*NN[n]*"_"*string(i)*".txt")[:,3])
    end
    h = StatsBase.fit(Histogram,log10.(PXA),nbins=200)
    x,y = [collect(h.edges[1]), h.weights]
    display(Plots.plot!(x[1:end-1],log10.(y) ./ x[1:end-1] ,label="N = "*string(NN[n])))
    # display(Plots.histogram!(log10.(PXA), normalize=:pdf,label=L"N = "*NN[n],yaxis=:log10))
end
xaxis!((-4,0))
xticks!([-4:1:0;], [L"10^{-4}",L"10^{-3}",L"10^{-2}",L"10^{-1}",L"10^0"])
yticks!([1], [""])
# yticks!([1:4;], reverse([L"10^{4}",L"10^{3}",L"10^{2}",L"10^{1}"]))
xlabel!(L"x_{\sigma}")
ylabel!(L"P(x_{\sigma})")
# Plots.pdf("convN")


# @ N=100 fixed, convM, Justification of the choice of M
    # N = 100
    # MM = [50,100,200,500,1000,1300,1600]
    # data = readdlm(".\\convM\\N100convM.txt")
    # data = reshape(data,(length(data),1))
    # pyplot()
    # Plots.plot()
    # for m in reverse(MM)
    #     PX = data[1:m*N]
        # h = StatsBase.fit(Histogram,log10.(PX),nbins=100)
        # x,y = [collect(h.edges[1]), h.weights]
        # display(Plots.plot!(x[1:end-1],log10.(y) .+ log10(m),label="M x N = "*string(m*N)))
    # end
    # xaxis!((-3,1))
    # xlabel!(L"x_σ")
    # ylabel!(L"P(x_σ)")
    # xticks!([-3:1:1;], [L"10^{-3}",L"10^{-2}",L"10^{-1}",L"10^0",L"10^1"])
    # yticks!([0:3;], reverse([L"10^{3}",L"10^{2}",L"10^{1}",L"10^0"]))
    # #Plots.pdf("convM")

# Plots Distribution Stopping Force
    # force_stoppA = zeros(M)
    # PXA = Float64[]
    # for i in 1:M
    #     global PXA = vcat(PXA,readdlm(".\\convN\\PX_k"*string(k_gamma)*"_MF"*string(k_MF)*"_N50_"*string(i)*".txt")[:,3])
    #     # force_stoppA[i] = readdlm(".\\convN\\stop_force"*string(k_gamma)*"_MF"*string(k_MF)*"_N50_"*string(i)*".txt")[1]
    # end
    # pyplot()
    # Plots.histogram(force_stoppA, label="Above", normalize=:pdf)
    # Plots.histogram!(force_stoppB, label="Below", normalize=:pdf)
    # xlabel!("f (critical force)")
    # ylabel!("P(f)")
    # Plots.plot(xA,yA,lw=5)
    # Plots.plot!(xB,yB,lw=5)
    # Plots.pdf("Distribution_stopping_force")


## Theta vs k_MF
    # kMF = [0.5,1,2,5,10,100]
    # thetaMF = [1.25,1.27,1.26,1.20,1.19,1.14]
    # yerr = 0.05
    # Plots.plot(kMF,thetaMF,xaxis=:log,yerr=yerr, marker=:circle, legend=nothing)
    # xlabel!(L"k_{MF}")
    # ylabel!("θ")
    # Plots.pdf("thetaTrendMF")

# ## Theta vs k_MF
    # kGamma = [0.5,1,1.67,2,3,6]
    # thetagamma = [1.19,1.27,1.24,1.24,1.19,1.25]
    # yerr = 0.05
    # Plots.plot(kGamma,thetagamma,yerr=yerr,marker=:circle, legend=nothing)
    # xlabel!(L"k_{Γ}")
    # ylabel!("θ")
    # Plots.pdf("thetaTrendGamma")
