using DelimitedFiles , Distributed , Plotly , PyPlot
println() ; println() ; println()

# Physical Parameters
    N = 100
    m   = 1 # mass of each particle (all blocks are identical for now)
    k_MF   = 1 # MF strength
    strength_potential = 1.0 # note that k_spring = 2 in our model
    #eta = 1/10*sqrt(4*m*(k_MF + 2.0*strength_potential)) # coeff 1/10 in front chosen so that one is sure to be in the underdamped regime
    eta = 0.34641016151377546
# Numerical Parameters
    dt = 0.1
    t_final = 5000
## Creation of the distribution
    distribution_name = "gamma"
    k_gamma = 1

# Param for distribution_param = 10, k=1
f_target_A  = vcat(0.0,10.0:0.2:23.0,25.0:5.0:80.0)
M_disorder_A = vcat(1,fill(20,78))
f_target_B = vcat(0.0,7.0:0.5:15.0,16.0:2.0:80.0)
M_disorder_B = vcat(1,fill(50,50))


force_A = rep_ID_force_A = rep_ID_dis_A = []
for i in 1:length(f_target_A)
    global force_A = vcat(force_A,fill(f_target_A[i],M_disorder_A[i]))
        # this creates [f1 ... f1 f2 ... f2 f3 ... f3 etc etc], each time, there is m_i times the force_A f_i
        # e.g : [1 2 2 2 4 5 5] if the associated numbers of disorder requiered are [1 3 1 2]
    global rep_ID_force_A = vcat(rep_ID_force_A,fill(string(f_target_A[i]),M_disorder_A[i]))
    global rep_ID_dis_A = vcat(rep_ID_dis_A,[string(s) for s in 1:M_disorder_A[i] ])
end
force_B = rep_ID_force_B = rep_ID_dis_B = []
for i in 1:length(f_target_B)
    global force_B = vcat(force_B,fill(f_target_B[i],M_disorder_B[i]))
        # this creates [f1 ... f1 f2 ... f2 f3 ... f3 etc etc], each time, there is m_i times the force_B f_i
        # e.g : [1 2 2 2 4 5 5] if the associated numbers of disorder requiered are [1 3 1 2]
    global rep_ID_force_B = vcat(rep_ID_force_B,fill(string(f_target_B[i]),M_disorder_B[i]))
    global rep_ID_dis_B = vcat(rep_ID_dis_B,[string(s) for s in 1:M_disorder_B[i] ])
end
#
lenf_A = length(force_A)
@time pmap(simulation_protocol_above,fill(NA,lenf_A),fill(dt,lenf_A),fill(t_final,lenf_A),fill(m,lenf_A),fill(eta,lenf_A),fill(k_MF,lenf_A),force_A,fill(distribution_name,lenf_A),fill(distribution_param,lenf_A),rep_ID_dis_A,rep_ID_force_A)

lenf_B = length(force_B)
@time pmap(simulation_protocol_below,fill(NB,lenf_B),fill(dt,lenf_B),fill(t_final,lenf_B),fill(m,lenf_B),fill(eta,lenf_B),fill(k_MF,lenf_B),force_B,fill(distribution_name,lenf_B),fill(distribution_param,lenf_B),rep_ID_dis_B,rep_ID_force_B)

path_A = ".\\protocol_above\\"
vel_SS_average_on_disorder_A = zeros(length(f_target_A))
std_SS_average_on_disorder_A = zeros(length(f_target_A))
ratio_average_on_disorder_A = zeros(length(f_target_A))
min_A = zeros(length(f_target_A))
max_A = zeros(length(f_target_A))
for i in 1:length(f_target_A)
    ID_force = string(f_target_A[i])
    vel_tmp_A = zeros(M_disorder_A[i])
    for config_number in 1:M_disorder_A[i]
        ID_config = string(config_number)

        ## Read velocity files
        fname_A = path_A*"MeanVelocityStep_"*ID_force*"_"*ID_config*".txt" ; vel_A = readdlm(fname_A)
        vel_tmp_A[config_number] = mean(vel_A[floor(Int,0.9*end):end])
    end
    min_A[i] = abs(minimum(vel_tmp_A)-mean(vel_tmp_A))
    max_A[i] = abs(maximum(vel_tmp_A)-mean(vel_tmp_A))
    vel_SS_average_on_disorder_A[i] = mean(vel_tmp_A)
    std_SS_average_on_disorder_A[i] = std(vel_tmp_A)
    ratio_average_on_disorder_A[i] = 100*length(vel_tmp_A[vel_tmp_A .< 0.1])/length(vel_tmp_A)  # percentage of stopped configurations
end

path_B = ".\\protocol_below\\"
vel_SS_average_on_disorder_B = zeros(length(f_target_B))
std_SS_average_on_disorder_B = zeros(length(f_target_B))
ratio_average_on_disorder_B = zeros(length(f_target_B))
min_B = zeros(length(f_target_B))
max_B = zeros(length(f_target_B))
for i in 1:length(f_target_B)
    ID_force = string(f_target_B[i])
    vel_tmp_B = zeros(M_disorder_B[i])
    for config_number in 1:M_disorder_B[i]
        ID_config = string(config_number)

        ## Read velocity files
        fname_B = path_B*"MeanVelocityStep_"*ID_force*"_"*ID_config*".txt" ; vel_B = readdlm(fname_B)
        vel_tmp_B[config_number] = mean(vel_B[floor(Int,0.9*end):end])
    end
    min_B[i] = abs(minimum(vel_tmp_B)-mean(vel_tmp_B))
    max_B[i] = abs(maximum(vel_tmp_B)-mean(vel_tmp_B))
    vel_SS_average_on_disorder_B[i] = mean(vel_tmp_B)
    std_SS_average_on_disorder_B[i] = std(vel_tmp_B)
    ratio_average_on_disorder_B[i] = 100*length(vel_tmp_B[vel_tmp_B .< 0.1])/length(vel_tmp_B)  # percentage of stopped configurations
    # pyplot()
    # display(Plots.histogram(vel_tmp_B,title=ID_force))
end




# # # Plots


pyplot()
Plots.plot(f_target_A,vel_SS_average_on_disorder_A,lw=2,box=true,label="Above",grid=false,ribbon=std_SS_average_on_disorder_A,fillalpha=.5)
Plots.plot!(f_target_B,vel_SS_average_on_disorder_B,yerr=(min_B,max_B),lw=2,box=true,label="Below",grid=false,ribbon=std_SS_average_on_disorder_B,fillalpha=.5)
Plots.plot!(1:10, NaN.*(1:10), label = "Stop Ratio Above [%]", linecolor=:blue, ls=:dash, grid=false, legend=:top)
Plots.plot!(1:10, NaN.*(1:10), label = "Stop Ratio Below [%]", linecolor=:red, ls=:dash ,grid=false)
Plots.plot!(Plots.twinx(),f_target_A,ratio_average_on_disorder_A,ls=:dash,lc=:blue,legend=false,ymirror=true,ylabel="ff")
Plots.plot!(Plots.twinx(),f_target_B,ratio_average_on_disorder_B,ls=:dash,lc=:red,legend=false)
xlabel!("Force")
ylabel!("Steady State Velocity")
Plots.pdf("Hysteresis_k1_continuous_protocol")

# fcA = 16.5
# fcB = 16.5
# indexA = length(f_target_A) - length(filter(x->x>fcA,f_target_A)) +1
# indexB = length(f_target_B) - length(filter(x->x>fcB,f_target_B)) +1
# Plots.plot(f_target_A[indexA:end] .- fcA,vel_SS_average_on_disorder_A[indexA:end],lw=2,label="Above", xaxis=:log, yaxis=:log,ribbon=std_SS_average_on_disorder_A[indexA:end],fillalpha=.5)
# Plots.plot!(f_target_B[indexB:end] .- fcB ,vel_SS_average_on_disorder_B[indexB:end],lw=2,label="Below", xaxis=:log, yaxis=:log,ribbon=std_SS_average_on_disorder_B[indexB:end],fillalpha=.5)
# # Plots.plot!([5:0.5:7.5], [14.2, 15.64,17.3,18.84,19.64,21.69],markershape=:circle,lw=2,label="CRH1")
# xlabel!("Force")
# ylabel!("Steady State Velocity")
# # Plots.pdf("Hysteresis_m01_loglog")
#
# # title!("NA = "*string(NA)*" & NB = "*string(NB))
# p2 = Plots.plot(f_target_A,ratio_average_on_disorder_A)
# Plots.plot(p1,p2,layout=(2,1))
#
if workers()!=[1]
    println("WARNING : ", length(workers()), " cores still in activity.")
end

# rmprocs(workers())
