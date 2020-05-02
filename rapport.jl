# ## Creation of the distribution
# distribution_name = "gamma"
# # config_disorder1 = config_disorder =  createCusps(distribution_name,1)
# # config_disorder05 =  createCusps(distribution_name,0.5)
# # config_disorder2 =  createCusps(distribution_name,2)
# # config_disorder5 =  createCusps(distribution_name,5)
#
# using Plots
# using PyPlot
# using DelimitedFiles
#
# using LaTeXStrings
# pyplot()
# # dx=0.1
# # x = [i for i in 1:dx:90]
# # force = [getForcePotential(element,config_disorder) for element in x]
# # Plots.plot(x,force,lw=2,label="Force Landscape",xlabel=L"$x$",ylabel=L"$f(x)$")
# # PyPlot.savefig("C:\\Users\\Ylann Rouzaire\\JupyterNotebook\\Internship Summer 2019\\images_report\\force.png")
#
# # pot = [getPotential(element,config_disorder) for element in x]
# # Plots.plot(x,pot,lw=2,label="Potential Landscape",xlabel=L"$x$",ylabel=L"$V(x)$")
# # PyPlot.savefig("C:\\Users\\Ylann Rouzaire\\JupyterNotebook\\Internship Summer 2019\\images_report\\potential.png")
#
# # Plots.histogram(config_disorder05,nbins=90,label="Shape = 0.5")
# # Plots.histogram!(config_disorder1,nbins=90,label="Shape = 1")
# # Plots.histogram!(config_disorder2,nbins=90,label="Shape = 2")
# # title!("Different shapes of the Gamma distribution")
# # Plots.savefig("C:\\Users\\Ylann Rouzaire\\JupyterNotebook\\Internship Summer 2019\\images_report\\hist.png")
#
# ## Convergence Tests
# # Physical Parameters
# N   = 1 # number of particles
# m   = 1 # mass of each particle (all blocks are identical for now)
# k_MF   = 1 # MF strength
# f = 0
# strength_potential = 1.0 # note that k_spring = 2 in our model
# eta = 1/10*sqrt(4*m*(k_MF + 2.0*strength_potential)) # coeff 1/10 in front chosen so that one is sure to be in the underdamped regime
# # Numerical Parameters
# dt = 0.1
# t_final = 50
#
# simulation_protocol_above() (N,dt,t_final,m,eta,k_MF,f,config_disorder,"gamma",1)
# fname = ".\\MeanVelocity.txt"; vel = readdlm(fname)
# fname = ".\\MeanPosition.txt"; pos = readdlm(fname)
#
# # Plots.plot(0:dt:t_final-dt,vel, label="f = 0")
# # xlabel!("Time")
# # ylabel!("v(t)")
# # Plots.savefig("C:\\Users\\Ylann Rouzaire\\JupyterNotebook\\Internship Summer 2019\\images_report\\velf0.png")
#
#
# Plots.plot(0:dt:t_final-dt,vel, label="f = 0")
# xlabel!("Time")
# ylabel!("v(t)")
# Plots.savefig("C:\\Users\\Ylann Rouzaire\\JupyterNotebook\\Internship Summer 2019\\images_report\\velf0.png")
#
#
# ## function of N
# f_target_A  = vcat(0, 3.1:0.1:5.3 , 5.5  ,12)
# f_target_B  = vcat(0, 5:0.1:7,7.5,8,12)
#
# NNA    = [ "N=50" "N=75" "N=100" "N=200" "N=300" "N=500"]
# NNB    = ["N=5" "N=10" "N=50" "N=100" "N=200" "N=300" "N=500"]
#
# vA5 = [2.8092e-5, 5.7829e-6, -0.00107178, -0.00202085, -0.000356033, 0.131192, 0.0632591, 0.943864, 2.93071, 2.46257, 1.56072, 3.44703, 5.30689, 7.37316, 10.2154, 8.53412, 10.3437, 11.8016, 12.8211, 13.0042, 13.7159, 13.864, 14.2073, 14.7803, 15.425, 33.5864]
# vA10 = [3.00153e-5, 0.000116068, 0.000328431, 0.0181187, 0.0102135, 0.166874, 0.773214, 0.757766, 1.93694, 3.06743, 4.7453, 6.10564, 6.4305, 8.18902, 9.44015, 10.4349, 12.771, 13.1543, 13.4429, 13.8024, 14.1629, 14.3866, 14.7286, 15.0044, 15.6395, 34.2002]
# vA50  = [2.44597e-5, -0.000149636, 0.000395859, 0.000564498, 0.0364846, 0.153567, 0.422223, 2.53163, 4.16127, 5.27925, 5.76465, 6.28734, 7.13778, 7.99697, 10.6103, 12.1438, 12.6746, 13.1768, 13.5456, 13.899, 14.2284, 14.5481, 14.8611, 15.1758, 15.7876, 34.6056]
# vA75 = [3.09693e-5, 0.000475321, 0.000614857, -0.000546867, 0.030234, 0.1301, 0.305951, 2.49271, 4.8754, 5.46518, 5.86179, 6.27649, 6.82264, 8.15265, 10.7747, 12.1828, 12.7644, 13.165, 13.5031, 13.8642, 14.2331, 14.5723, 14.8663, 15.1793, 15.785, 34.6487]
# vA100 = [1.84012e-5, -2.7029e-5, -0.000191335, 0.000465035, 0.0298615, 0.0381255, 0.836168, 2.78613, 4.82793, 5.40786, 5.84362, 6.2193, 6.60644, 7.66064, 10.1713, 12.1448, 12.8071, 13.2347, 13.5597, 13.9426, 14.2724, 14.5787, 14.8956, 15.1841, 15.7738, 34.6359]
# vA200 = [2.34712e-5, -0.000286519, 0.000369819, -0.000467173, 0.00530396, 0.0709413, 0.304617, 2.68039, 4.97276, 5.49124, 5.83396, 6.21919, 6.8948, 8.06807, 10.9087, 12.2649, 12.8241, 13.2105, 13.5613, 13.9155, 14.2434, 14.5711, 14.9, 15.1854, 15.8052, 34.6743]
# vA300 = [2.59421e-5, 4.90894e-5, 0.000382743, 0.000372684, 0.0141703, 0.140376, 0.719432, 3.03345, 4.99259, 5.48788, 5.85332, 6.2971, 6.88019, 7.83549, 10.8092, 12.2252, 12.7438, 13.2132, 13.5879, 13.9276, 14.2635, 14.5797, 14.8931, 15.1951, 15.8061, 34.6746]
# vA500 = [3.14036e-5, -0.000291815, 0.000273429, -0.000251656, 0.00736, 0.107666, 0.628853, 3.50668, 5.00694, 5.47914, 5.86035, 6.22661, 6.8042, 7.89868, 11.0246, 12.3105, 12.8094, 13.2355, 13.5904, 13.9288, 14.2606, 14.5786, 14.9, 15.1989, 15.8238, 34.6873]
#
# vB5 = [-5.70588e-14, -9.06015e-5, -7.83062e-5, -8.02166e-6, -5.75013e-6, 1.16365e-6, -1.68378e-6, -6.61763e-6, -1.40007e-5, -1.19601e-6, -8.57677e-5, -7.835e-5, -5.2788e-6, 1.74383, -5.33655e-5, -5.16599e-5, -5.157e-5, 1.84862, 1.88074, 3.83417, 3.875, 3.93523, 21.1351, 22.774, 33.6584]
# vB10 = [7.34304e-15, -0.000104418, -0.000106483, -0.000108545, -0.000110606, -8.16606e-5, -8.2716e-5, -2.02012e-6, -3.01553e-5, -6.00492e-5, -6.65216e-5, -7.41968e-5, -2.41767e-5, -1.61278e-5, 1.76986, 2.16399e-5, 3.41393e-5, 1.87404, 1.89937, 3.86557, 3.94436, 11.949, 21.346, 22.7833, 34.1978]
# vB50  = [1.93382e-14, -4.2363e-5, -2.52678e-5, -1.99375e-5, -1.90088e-5, -1.45307e-5, -1.33336e-5, 1.70378e-5, 1.47078, 1.62366, 1.95075e-5, 1.71497, 1.74753, 3.5519, 7.21003, 9.1851, 9.34886, 13.2632, 15.3774, 17.59, 17.8657, 18.1197, 21.5348, 23.0976, 34.5525]
# vB100 = [-8.8345e-15, -5.69381e-5, -4.98472e-5, -5.42586e-5, -3.396e-5, -1.51692e-5, 1.37209, 3.03686, 3.13349, 3.27148, 6.71958, 8.45637, 10.3748, 15.9334, 18.0332, 18.402, 18.7085, 19.0091, 19.2976, 19.5961, 19.9015, 20.1847, 21.593, 23.1044, 34.6349]
# vB200 = [3.80321e-17, -5.51613e-5, -5.30683e-5, -4.67104e-5, -3.56419e-5, -4.88788e-5, -3.24137e-5, -3.28323e-5, -7.64078e-5, -1.56621e-5, 1.48796, 6.50101, 8.64673, 12.1374, 14.3595, 16.5228, 16.8332, 17.0983, 19.288, 19.601, 19.8949, 20.1894, 21.637, 23.0934, 34.6823]
# vB300 = [4.02374e-15, -5.96943e-5, -4.99385e-5, -4.88048e-5, -4.1533e-5, -4.88125e-5, -4.14041e-5, -2.65068e-5, -1.09553e-5, 0.000710036, 3.20248, 5.08202, 8.51733, 13.8537, 18.0264, 18.4072, 18.7251, 19.0271, 19.3221, 19.6103, 19.9109, 20.1947, 21.668, 23.1137, 34.6682]
# vB500 = [1.68494e-15, -4.75093e-5, -4.43956e-5, -3.91288e-5, -3.23601e-5, -3.52069e-5, -2.07514e-5, -1.51665e-5, -5.14753e-6, 1.22818e-5, 2.712, 9.82841,13.8475, 15.8801, 17.816, 18.411, 18.7235, 19.0239, 19.3223, 19.6134, 19.9125, 20.2033, 21.6506, 23.1115, 34.657]
#
# Plots.plot(f_target_A,[vA50,vA75,vA100,vA200,vA300,vA500],lw=2, label=NNA,title="Above")
# xlabel!("Force")
# Plots.ylabel!(L"overline{v_{SS}}")
# Plots.plot(f_target_B,[vB5,vB10,vB50,vB100,vB200,vB300,vB500],lw=2, label=NNB,title="Below")


## Figures force velocity step Continuous

# # Physical Parameters
#     N   = 100 # number of particles
#     mass  = 1 # mass of each particle (all blocks are identical for now)
#     k_MF   = 1 # MF strength
#     strength_potential = 1.0 # note that k_spring = 2 in our model
#     #eta = 1/10*sqrt(4*m*(k_MF + 2.0*strength_potential)) # coeff 1/10 in front chosen so that one is sure to be in the underdamped regime
#     eta = 0.34641016151377546
# # Numerical Parameters
#     dt = 0.1
#     t_final = 1500
# ## Creation of the distribution
#     distribution_name = "gamma"
#     distribution_param = 10
#     target_force = 20
#
# simulation_protocol_above(N,dt,t_final,mass,eta,k_MF,target_force,distribution_name,distribution_param)
# simulation_protocol_below(N,dt,t_final,mass,eta,k_MF,target_force,distribution_name,distribution_param)
# protocol_continuous_above(N,dt,t_final,mass,eta,k_MF,target_force,distribution_name,distribution_param)
# protocol_continuous_below(N,dt,t_final,mass,eta,k_MF,target_force,distribution_name,distribution_param)
#
# fname_B = ".\\protocol_above\\Force_step.txt"; stepA = readdlm(fname_B)
# fname_B = ".\\protocol_below\\Force_step.txt"; stepB = readdlm(fname_B)
# fname_B = ".\\protocol_below\\force.txt"; contB = readdlm(fname_B)
# fname_B = ".\\protocol_above\\force.txt"; contA = readdlm(fname_B)
#
#
# fname_B = ".\\protocol_above\\MeanVelocity__.txt"; velA = readdlm(fname_B)
# fname_B = ".\\protocol_below\\MeanVelocity__.txt"; velB = readdlm(fname_B)
# fname_B = ".\\protocol_below\\MeanVelocityContinuous__.txt"; velcB = readdlm(fname_B)
# fname_B = ".\\protocol_above\\MeanVelocityContinuous__.txt"; velcA = readdlm(fname_B)
#
#
# pyplot()
# Plots.plot(0:0.1:t_final,stepA,label="Step Protocol",lw=2)
# Plots.plot!(0:0.1:t_final,contA,label="Continuous Protocol",lw=2)
# Plots.xlabel!("t")
# Plots.ylabel!("f(t)")
# Plots.pdf("Force_protocols_Above")
#
# Plots.plot(0:0.1:t_final,stepB,label="Step Protocol",lw=2)
# Plots.plot!(0:0.1:t_final,contB,label="Continuous Protocol",lw=2)
# Plots.xlabel!("t")
# Plots.ylabel!("f(t)")
# Plots.pdf("Force_protocols_Below")
#
# Plots.plot(0:0.1:t_final,velB,label="Step Protocol",lw=2)
# Plots.plot!(0:0.1:t_final,velcB,label="Continuous Protocol",lw=2)
# Plots.xlabel!("t")
# Plots.ylabel!("v(t)")
# Plots.pdf("Velocity_protocols_Below")
#
# Plots.plot(0:0.1:t_final,velA,label="Step Protocol",lw=2)
# Plots.plot!(0:0.1:t_final,velcA,label="Continuous Protocol",lw=2)
# Plots.xlabel!("t")
# Plots.ylabel!("v(t)")
# Plots.pdf("Velocity_protocols_Above")
