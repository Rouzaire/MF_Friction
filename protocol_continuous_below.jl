@everywhere function protocol_continuous_below(N,dt,t_final,m,eta,k_MF,target_force,distrib,ID_config="",ID_force="")
let

t = 0.0
particle_position = fill(1e-4,N)
particle_velocity = zeros(N)

step_time = [50,100,t_final*0.75,t_final]

f_vel = open(".\\protocol_below\\MeanVelocityContinuous_"*ID_force*"_"*ID_config*".txt", "w")
f_force = open(".\\protocol_below\\force.txt", "w")

# Phase 1 : impose large force so that the system looses memory of initial conditions
force = 5 # arbitrary
while t < step_time[1]
    t = t + dt
    write(f_force,string(force) * "\n")
    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,m,force,eta,k_MF,distrib[j],j)
    end
end

# Phase 2 : stop system and wait 50s for the transiens to fade out
force = 0
while t < step_time[2]
    t = t + dt
    write(f_force,string(force) * "\n")
    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,m,force,eta,k_MF,distrib[j],j)
    end
end

# Phase 3 : increase force linearly at each timestep (approximation of continuity)
df = target_force*dt/(step_time[3] - step_time[2])
while t < step_time[3]
    force = force + df
    t = t + dt
    write(f_force,string(force) * "\n")
    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,m,force,eta,k_MF,distrib[j],j)
    end
end

# Phase 4 : wait until t_final with force constant (f=target_force by construction)
while t < step_time[4]
    t = t + dt
    write(f_force,string(force) * "\n")
    write(f_vel,string(getMeanVelocity(particle_velocity)) * "\n")
    for j in 1:N
        particle_position[j],particle_velocity[j] = VV(particle_position,particle_velocity[j],dt,m,force,eta,k_MF,distrib[j],j)
    end
end
close(f_vel)
close(f_force)

end # let
end # function
