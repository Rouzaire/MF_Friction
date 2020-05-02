# vel1 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_1.txt")
# vel2 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_2.txt")
# vel3 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_3.txt")
# vel4 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_4.txt")
# vel5 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_5.txt")
# vel6 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_6.txt")
# vel7 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_7.txt")
# vel8 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_8.txt")
# vel9 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_9.txt")
# vel10 = readdlm(".\\protocol_above\\MeanVelocityContinuous_8.5_10.txt")
#
# pyplot()
# Plots.plot([0:dt:t_final],vel1)
# Plots.plot!([0:dt:t_final],vel2)
# Plots.plot!([0:dt:t_final],vel3)
# Plots.plot!([0:dt:t_final],vel4)
# Plots.plot!([0:dt:t_final],vel5)
# Plots.plot!([0:dt:t_final],vel6)
# Plots.plot!([0:dt:t_final],vel7)
# Plots.plot!([0:dt:t_final],vel8)
# Plots.plot!([0:dt:t_final],vel9)
# Plots.plot!([0:dt:t_final],vel10)
#
#
# vel1 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_1.txt")
# vel2 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_2.txt")
# vel3 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_3.txt")
# vel4 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_4.txt")
# vel5 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_5.txt")
# vel6 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_6.txt")
# vel7 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_7.txt")
# vel8 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_8.txt")
# vel9 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_9.txt")
# vel10 = readdlm(".\\protocol_below\\MeanVelocityContinuous_11.0_10.txt")
#
# pyplot()
# Plots.plot([0:dt:t_final],vel1)
# Plots.plot!([0:dt:t_final],vel2)
# Plots.plot!([0:dt:t_final],vel3)
# Plots.plot!([0:dt:t_final],vel4)
# Plots.plot!([0:dt:t_final],vel5)
# Plots.plot!([0:dt:t_final],vel6)
# Plots.plot!([0:dt:t_final],vel7)
# Plots.plot!([0:dt:t_final],vel8)
# Plots.plot!([0:dt:t_final],vel9)
# Plots.plot!([0:dt:t_final],vel10)
# global t = 0.0
# dt = 0.1
# N = 10
# distrib = [createCusps("gamma",1) for n in 1:N]
# particle_position = fill(1e-4,N)
# particle_velocity = zeros(N)
# noise_all_particles = zeros(N)
# mass = force = k_MF = 1
# eta = 0.35
# for tt in 1:10
#     for j in 1:N
#         particle_position[j],particle_velocity[j],noise_all_particles[j] = VV(particle_position,particle_velocity,noise_all_particles[j],dt,mass,force,eta,k_MF,distrib[j],j)
#     end
# end
