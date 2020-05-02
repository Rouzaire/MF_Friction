@everywhere using DelimitedFiles , Distributed , Plots, LaTeXStrings, StatsPlots, StatsBase
# Physical Parameters
    N   = 100 # number of particles
    mass   = 1 # mass of each particle (all blocks are identical for now)
    k_MF   = 1 # MF strength
    eta = 0.34641016151377546
    force = 1
# Numerical Parameters
    dt = .1
    t_final = 10
    t = 0.0
    particle_position = fill(1e-3,N)
    particle_velocity = zeros(N)
## Creation of the distribution
    distribution_name = "gamma"
    distribution_param = 1 # k_gamma
    distrib = [createCusps(distribution_name,distribution_param) for n in 1:N]
## Output files
f_pos = open(".\\Basic Simulation\\MeanPosition.txt", "w")
f_vel = open(".\\Basic Simulation\\MeanVelocity.txt", "w")
f_std = open(".\\Basic Simulation\\std.txt", "w")
f_eca = open(".\\Basic Simulation\\ecartmax.txt", "w")
f_acc = open(".\\Basic Simulation\\acc.txt", "w")
f_damp = open(".\\Basic Simulation\\damping.txt", "w")
f_kMF = open(".\\Basic Simulation\\kMF.txt", "w")
f_pot = open(".\\Basic Simulation\\force_potential.txt", "w")
f_meanacc = open(".\\Basic Simulation\\meanacc.txt", "w")

## Simulation
while t < t_final
    global t = t + dt
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
    end
    write(f_pos,string(getMeanPosition(particle_position)),"\n")
    write(f_vel,string(getMeanVelocity(particle_velocity)),"\n")
    write(f_std,string(getStd(particle_position)),"\n")
    write(f_eca,string(maximum(particle_position) - minimum(particle_position)),"\n")
    write(f_meanacc,string(getMeanAcceleration(particle_position,particle_velocity,mass,force,eta,k_MF,distrib)),"\n")
    write(f_acc,string(acceleration(particle_position[1],particle_velocity[1],mass,force,eta,k_MF,distrib[1],getMeanPosition(particle_position))),"\n")

    write(f_damp,string(getForceDamping(particle_velocity[1],eta)),"\n")
    write(f_pot,string(getForcePotential(particle_position[1],distrib[1])),"\n")
    write(f_kMF,string(getForceMF(particle_position,1,k_MF)),"\n")

end
close(f_pos)
close(f_vel)
close(f_std)
close(f_eca)
close(f_meanacc)
close(f_acc)
close(f_damp)
close(f_kMF)
close(f_pot)


## Data Analysis
mean_pos = readdlm(".\\Basic Simulation\\MeanPosition.txt")
mean_vel = readdlm(".\\Basic Simulation\\MeanVelocity.txt")
stddev = readdlm(".\\Basic Simulation\\std.txt")
ecartmax = readdlm(".\\Basic Simulation\\ecartmax.txt")
meanacc = readdlm(".\\Basic Simulation\\meanacc.txt")
acc1 = readdlm(".\\Basic Simulation\\acc.txt")

feta = readdlm(".\\Basic Simulation\\damping.txt")
fMF = readdlm(".\\Basic Simulation\\kMF.txt")
fpot = readdlm(".\\Basic Simulation\\force_potential.txt")

# p0 = mean_pos
# v0 = mean_vel
pyplot()
# # Plots.plot!([0:dt:t_final-dt],mean_pos,label="f=5")
# Plots.plot!([0:dt:t_final-dt],mean_vel,label="f=1")
# xaxis!(L"t")
# # yaxis!(L"\overline{x}(t)")
# yaxis!(L"\overline{v}(t)")
# Plots.pdf("vbarN100")
# Plots.pdf("vN1")

#
#
#
# Plots.histogram!(particle_position/getMeanPosition(particle_position),normalize=:pdf,label="f=0")
# xaxis!(L"\dfrac{x}{\ \overline{x}\ } \ (t_{final})",(0,2))
# yaxis!(L"P(\dfrac{x}{\ \overline{x} \ } \ (t_{final}))")
# Plots.pdf("Distribution_positionsN100")

# Plots.histogram!(particle_velocity,normalize=:pdf,label="f=1")
# xaxis!(L"v(t_{final})")
# yaxis!(L"P(v(t_{final}))")
# Plots.pdf("Distribution_velocitiesN100")
# Plots.histogram(particle_velocity/getMeanVelocity(particle_velocity),nbins=20)
# Plots.plot([0:dt:t_final],-eta*mean_vel)
# t0 = stddev
# s0 = stddev./mean_pos
# t1 = stddev
# s1 = stddev./mean_pos
# t2 = stddev
# s2 = stddev./mean_pos
# t3 = stddev
# s3 = stddev./mean_pos
# t4 = stddev
# s4 = stddev./mean_pos
# Plots.plot([0:dt:t_final-dt],t0,yaxis=:log,label=L"k_{MF} = 0")
# Plots.plot!([0:dt:t_final-dt],t2,yaxis=:log,label=L"k_{MF} = 1")
# Plots.plot!([0:dt:t_final-dt],t3,yaxis=:log,label=L"k_{MF} = 10")
# Plots.plot!([0:dt:t_final-dt],t4,yaxis=:log,label=L"k_{MF} = 10^2")
# xlabel!("t")
# ylabel!(L"$σ(t)\, /\, \overline{x}(t)$")
# ylabel!(L"$σ(t)$")
# # title!("Impact of MF strength on Standard Deviation.")
# Plots.pdf("Std_over_time")
# #
# Plots.plot([i for i in 0:dt:t_final-dt],fill(force,Int(t_final/dt +1)),label=L"f=0")
# Plots.plot!([i for i in 0:dt:t_final-dt],feta,label=L"-\eta\,\dot{x}")
# Plots.plot!([i for i in 0:dt:t_final-dt],fMF,label=L"k_{MF} (\bar{x} - x)")
# Plots.plot!([i for i in 0:dt:t_final-dt],fpot,label=L"-V'(x)")
# xlabel!("t")
# ylabel!("Force")
# # title!("Force Budget over the trajectory")
# Plots.pdf("Force_Budget_Move")
