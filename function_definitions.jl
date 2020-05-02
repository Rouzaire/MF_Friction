using Distributed
@everywhere using Distributions, DelimitedFiles, Plots, Statistics, DSP

@everywhere function createCusps(Distribution_name,argument) # returns an Array of the location of the cusps, according to the chosen distribution
    if Distribution_name == "gamma" # Gamma(shape, scale)=Gamma(k, theta)
        # in this case, the argument must correspond to the k (shape factor)
        mean = 1 # arbitrary choice, constraining the second argument theta (scale factor) to be equal to mean/k. Note that the variance therefore becomes mean^2/k
        size = floor(Int,1e5)
        # has to be large enough so that the simulation stops before it reaches the end
            # but its length is one of the most computation

        epsy = 0.001 .+ rand(Gamma(argument,mean/argument),size) # 0.001 shift so that a 0-width well is not possible.
        epsy = vcat(0.0,cumsum(epsy))
        return epsy

    elseif Distribution_name == "crystal"
        mean = 1 # arbitrary choice
        size = 1e5 # has to be large enough so that the simulation stops before it reaches the end
        phase = rand(mean)
        epsy =  fill(mean,Int(size))
        epsy = vcat(0.0,cumsum(epsy)) .- phase
        return epsy

    elseif Distribution_name == "onewell"
        return [0,argument]

    elseif Distribution_name == "stopwell"
        tmp = cumsum(fill(10.0,argument))
        return vcat(0,tmp,tmp[end]+50)
    end
    end

@everywhere function getForcePotential(position,potential_landscape_1particle) # returns -V'(eps) = -2(eps-Delta) at a given position
    cusp_before_index = searchsortedfirst(potential_landscape_1particle,position) - 1  # should be the complexity-optimal build-in function
    cusp_after_index = cusp_before_index + 1
    return (-2.0*(position - 0.5*(potential_landscape_1particle[cusp_after_index] + potential_landscape_1particle[cusp_before_index])))
            # since eps = position - eps_before and 2*Delta = eps_after - eps_before
end

@everywhere function getForceMF(position_all_particles,particleID,k_MF) # returns the force due to the Mean Field Coupling : k_MF(x_bar - x)
        return (k_MF*(getMeanPosition(position_all_particles) - position_all_particles[particleID]))
end

@everywhere function getForceDamping(velocity_1particle,eta)  # returns the force due to the damping : -eta*velocity
        return (-eta*velocity_1particle)
end

@everywhere function getWellNumber(position,potential_landscape_1particle) # returns the number of the well at a given position
    return (searchsortedfirst(potential_landscape_1particle,position) - 1)
end

@everywhere function getDelta(position,potential_landscape_1particle) # returns the half width (hence the 0.5 factor) Delta of the well at a given position
    cusp_before_index = searchsortedfirst(potential_landscape_1particle,position) - 1
    cusp_after_index = cusp_before_index + 1
    return 0.5*(potential_landscape_1particle[cusp_after_index] - potential_landscape_1particle[cusp_before_index])
end

@everywhere function getDelta0PlusMinus(position,potential_landscape_1particle) # for a given position, returns a 3x1 Array with the half widths of the well before (DeltaMinus), of the current well (Delta0), and of the next well (DeltaPlus)
    cusp_before_index = searchsortedfirst(potential_landscape_1particle,position) - 1  # should be the complexity-optimal build-in function
    cusp_after_index = cusp_before_index + 1
    DeltaMinus = 0.5*(potential_landscape_1particle[cusp_after_index - 1] - potential_landscape_1particle[cusp_before_index - 1])
    Delta0 = 0.5*(potential_landscape_1particle[cusp_after_index] - potential_landscape_1particle[cusp_before_index])
    DeltaPlus = 0.5*(potential_landscape_1particle[cusp_after_index + 1] - potential_landscape_1particle[cusp_before_index + 1])
    return [DeltaMinus,Delta0,DeltaPlus]
end

@everywhere function getPotential(position,potential_landscape_1particle) # returns V(eps) = (eps-Delta)^2 - Delta^2 at a given position eps
    cusp_before_index = searchsortedfirst(potential_landscape_1particle,position) - 1
    cusp_after_index = cusp_before_index + 1
    delta = 0.5*(potential_landscape_1particle[cusp_after_index] - potential_landscape_1particle[cusp_before_index])
    return ( (position - potential_landscape_1particle[cusp_before_index] - delta)^2 - delta^2 )
end


@everywhere function getForceYield(position,potential_landscape_1particle) # returns the force to yielding = 4 Delta - 2*(position)
    eps_before_index = searchsortedfirst(potential_landscape_1particle,position) -1
    eps_after_index = eps_before_index + 1
    delta = 0.5*(potential_landscape_1particle[eps_after_index] - potential_landscape_1particle[eps_before_index])
    return ( 2*delta + getForcePotential(position,potential_landscape_1particle))
end


@everywhere function getMeanPosition(position_all_particles)
    return mean(position_all_particles)
end

@everywhere function getMeanVelocity(velocity_all_particles)
    return mean(velocity_all_particles)
end

@everywhere function getMeanAcceleration(position_all_particles, velocity_all_particles, mass,f,eta,k,potential_landscape_all_particles)
    mean_position = mean(position_all_particles)
    acc = [acceleration(position_all_particles[n],velocity_all_particles[n],mass,f,eta,k,potential_landscape_all_particles[n],mean_position) for n in 1:length(position_all_particles)]
    return mean(acc)
end

@everywhere function getStd(position_all_particles) # returns the standard deviation RESCALED BY THE MEAN POSITION
    return Statistics.std(position_all_particles)/getMeanPosition(position_all_particles)
end

@everywhere function acceleration(position_1particle,velocity_1particle,mass,f,eta,kMF,potential_landscape_1particle,mean_position) # without additionnal noise
    g = 1 # elastic modulus
    return (-eta*velocity_1particle + f + g*getForcePotential(position_1particle,potential_landscape_1particle) + kMF*(mean_position - position_1particle))/mass
end

@everywhere function acceleration(position_1particle,velocity_1particle,mass,f,eta,kMF,potential_landscape_1particle,mean_position,noise) # WITH additionnal noise
    g = 1 # elastic modulus
    return (-eta*velocity_1particle + f + g*getForcePotential(position_1particle,potential_landscape_1particle) + kMF*( (mean_position - position_1particle) + noise ) )/mass
end

@everywhere function VV(position_all_particles,velocity_1particle,dt,mass,f,eta,k,potential_landscape_1particle,particleID) # Velocity Verlet without additionnal noise
    mean_position = Statistics.mean(position_all_particles)
    x = position_all_particles[particleID]
    a_before = acceleration(x,velocity_1particle,mass,f,eta,k,potential_landscape_1particle,mean_position)
    position_after = x + dt*velocity_1particle + 0.5*dt*dt*a_before
    velocity_estimate = velocity_1particle + 0.5*dt*(a_before + acceleration(position_after,velocity_1particle + dt*a_before,mass,f,eta,k,potential_landscape_1particle,mean_position))
        # one needs this velocity estimate because the a(t+dt) depends on v(t+dt) (cyclic problem).
        # In a first time, one estimates v(t+dt) roughly with Euler
    velocity_after = velocity_1particle + 0.5*dt*(a_before + acceleration(position_after,velocity_estimate,mass,f,eta,k,potential_landscape_1particle,mean_position))

    return [position_after,velocity_after]
end

@everywhere function VV(position_all_particles,velocity_all_particles,noise_1particle,dt,mass,f,eta,k,potential_landscape_1particle,particleID) # Velocity Verlet WITH additionnal noise
    mean_position = Statistics.mean(position_all_particles)
    mean_velocity = Statistics.mean(velocity_all_particles)

    # Creation of the noise
    tau = sqrt(2/mass) # 2 = 2*elastic_modulus = 2*1. It is the value of k in sqrt(k/m) = the typical timescale in a damped HO
    c0 = 0 # coeff to tune the variance of the following Gaussian
    noise = noise_1particle*(1-dt/tau) + rand(Normal(0,c0*mean_position*2/tau/dt))

    x = position_all_particles[particleID]
    v = velocity_all_particles[particleID]
    a_before = acceleration(x,v,mass,f,eta,k,potential_landscape_1particle,mean_position,noise)
    position_after = x + dt*v + 0.5*dt*dt*a_before
    velocity_estimate = v + 0.5*dt*(a_before + acceleration(position_after,v + dt*a_before,mass,f,eta,k,potential_landscape_1particle,mean_position,noise))
        # one needs this velocity estimate because the a(t+dt) depends on v(t+dt) (cyclic problem).
        # In a first time, one estimates v(t+dt) roughly with Euler
    velocity_after = v + 0.5*dt*(a_before + acceleration(position_after,velocity_estimate,mass,f,eta,k,potential_landscape_1particle,mean_position,noise))

    return [position_after,velocity_after,noise]
end
