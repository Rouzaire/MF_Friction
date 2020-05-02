@everywhere function protocol_continuous_above(N,dt,t_final,mass,eta,k_MF,target_force,distrib,ID_config="",ID_force="")
let

t = 0.0
particle_position = fill(1e-4,N)
particle_velocity = zeros(N)

f_initial = 2*target_force
step_time = [50,t_final*0.75,t_final]

f_vel = open(".\\protocol_above\\MeanVelocityContinuous_"*ID_force*"_"*ID_config*".txt", "w")
f_force = open(".\\protocol_above\\force.txt", "w")
f_std = open(".\\protocol_above\\std.txt", "w")
f_acc = open(".\\protocol_above\\MeanAcc.txt", "w")

# Phase 1 : wait for transients to fade out at f = f_initial
force = f_initial
while t < step_time[1]
    t = t + dt
    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    write(f_force,string(force) * "\n")
    write(f_std,string(getStd(particle_position)) * "\n")
    # write(f_acc,string(mean([acceleration(x,v,m,f,eta,k,epsy_1particle)])) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
    end
end

# Phase 2 : decrease force linearly at each timestep (approximation of continuity)
df = (f_initial - target_force)*dt/(step_time[2] - step_time[1])
while t < step_time[2]
    force = force - df
    t = t + dt
    write(f_force,string(force) * "\n")

    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
    end
end

# Phase 3 : wait until t_final with force constant (f=target_force by construction)
while t < step_time[3]
    t = t + dt
    # write(f_vel,string(force) * "\n")
    write(f_force,string(force) * "\n")

    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,mass,force,eta,k_MF,distrib[j],j)
    end
end

close(f_vel)
close(f_force)
close(f_std)
close(f_acc)
end # let
end # function
