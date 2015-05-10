# This model refers to the paper: {place holder for the reference}.
# Interested readers are referred to the paper for further description of the model.

# Define packages to be used:
using JuMP # Optimization package
using JuMPChance # Chance constraints package
using Gurobi # Solver, could be altered if needed
using MAT # Package to interface with .mat files

# Include the file with input data functions
include("input.jl")
# Include the file recalculating the B matrix
include("matrix.jl")

# Main function solving chance constrained opf
function solve_ccopf(logfilename, refbus, line_probability_threshold, gen_probability_threshold, linebufferamt, mvaBase, loadscale, thermalLimitscale, buses, lines, farms, generatorlist, Binvb, Bhatsp, sumvar, hydro_idx, nuclear_idx, hydro_limit, ramp, robust_budget, loaded_lines_idx, loaded_lines_probability_threshold)
    # Dimensions of sets
    numbuses = length(buses) # Number of buses
    numlines = length(lines) # Number of lines
    numfarms = length(farms) # Number of farms

    const Γ=iround(robust_budget*numfarms) # Budget of uncertainty
    const VoLL = 10_000 # Value of lost load
    
    # Define the name of the model and solver settings
    m = ChanceModel(solver=GurobiSolver(Method=1,BarHomogeneous=1))

    # Define variables: 
    # Non-negative hourly average power output bounded by the maximum limit on each generator:
    @defVar(m, 0 <= pbar[i=generatorlist] <= buses[i].Pgmax) 
    # Regulation participation factor between 0 and 1:
    @defVar(m, 0 <= alpha[i=generatorlist] <= ((buses[i].Pgmax == 0) ? 0 : 1)) 
    # Power flow in each transmission line should be within the limits:
    @defVar(m, -lines[i].u <= barf[i=1:numlines] <= lines[i].u)

    @defVar(m, beta[1:(numbuses-1)])
    @defVar(m, delta[1:(numbuses-1)])
    @defVar(m, g[1:(numbuses-1)]) # generator variables
    @defVar(m, h[1:(numbuses-1)]) # \hat B^-1 \bar p
    @defVar(m, θ[1:(numbuses-1)]) # thetahat

    @defVar(m, slack_bus[1:(numbuses-1)] >= 0)

    # Define random variable for each wind farm with a given mean and variance
    @defIndepNormal(m, winderr[i=1:numfarms], mean=farms[i].deviation_mean_interval, var=farms[i].deviation_sigmasq_interval)

    # Set up the reference bus:
    @assert refbus == numbuses
    @assert !(refbus in generatorlist)

    # Define constraints on participation factors:
    for i in 2:length(hydro_idx)
        # Enforce uniform participation factors on hydro power plants
        @addConstraint(m, alpha[hydro_idx[i-1]] == alpha[hydro_idx[i]]) 
    end
    # Set participation factors of nuclear generators to 0 (provide no regulation due to the ramping/cycling limitations)
    @addConstraint(m, nuclear_alpha_fixed[i=nuclear_idx], alpha[i] == 0)
    # Integrality constraint on participation factors:
    @addConstraint(m, sumalpha, sum(alpha) == 1)

    # Exogenous outputs of hydro and nuclear generators:
    # Enforce the hourly average power output (exogenous value, obtained from a separate hydro scheduling procedure) of wind farms 
    @addConstraint(m, hydro_ub[k=1:length(hydro_idx)], pbar[hydro_idx[k]] <= hydro_limit[k])
    # Enforce the hourly average power output of nuclear power plants (always committed and output at their maximum power output - "must run" units) 
    @addConstraint(m, nuclear_fixed[i=nuclear_idx], pbar[i] == buses[i].Pgmax)

    # System-wide power balance constraint:    
    sumload = sum([getLoad(b) for b in buses]) # system-wide load
    sumloadminusmeanwind = sum([getLoadMinusMeanWind(b) for b in buses]) # system-wide net load, i.e. the system-wide load minus the system-wide wind
    # Constraint on power balance:
    @addConstraint(m, balance, sum(pbar) == sumloadminusmeanwind)

    Brow = Bhatsp'

    # \hat B delta = \alpha
    @addConstraint(m, defalpha[i=1:(numbuses-1)],
        sum{Brow.nzval[idx]*delta[Brow.rowval[idx]], idx in Brow.colptr[i]:(Brow.colptr[i+1]-1)} -
        (isgen(buses[i]) ? alpha[i] : 0) == 0)

    # \hat B h = p
    @addConstraint(m, defh[i=1:(numbuses-1)],
        sum{Brow.nzval[idx]*h[Brow.rowval[idx]], idx in Brow.colptr[i]:(Brow.colptr[i+1]-1)} -
        (isgen(buses[i]) ? pbar[i] : 0) - slack_bus[i] == 0)

    # theta = \hat B^{-1}(w - d) + h = \hat B^{-1}(p + w - d)
    @addConstraint(m, BB[i=1:(numbuses-1)], θ[i]-h[i] == Binvb[i])

    # Calculating flows in transmission lines:
    @addConstraint(m, flowangle[i=1:numlines], barf[i] +
        (lines[i].tail != numbuses ? -lines[i].y*θ[lines[i].tail] : 0) +
        (lines[i].head != numbuses ? lines[i].y*θ[lines[i].head] : 0) == 0)


   # Chance constraints on power flow limits
    for i in 1:numlines
        # Calculate expr to be used in chance constraints
        ccexpr = 0
        tail = lines[i].tail
        head = lines[i].head
        for j in 1:numfarms
            fexpr = 0
            if tail != numbuses
                fexpr -= delta[tail] - farms[j].colbarinv[tail]
            end
            if head != numbuses
                fexpr += delta[head] - farms[j].colbarinv[head]
            end
            ccexpr += fexpr*winderr[j]
        end
        prob_threshold = (i in loaded_lines_idx) ? loaded_lines_probability_threshold : line_probability_threshold

        # Formulate chance constraints for the + and - transmission limits:
        addConstraint(m, barf[i] + lines[i].y*ccexpr >= getThermalCapacity(lines[i], mvaBase)*(1-linebufferamt), with_probability=prob_threshold, uncertainty_budget_mean=Γ, uncertainty_budget_variance=Γ)
        addConstraint(m, barf[i] + lines[i].y*ccexpr <= -getThermalCapacity(lines[i], mvaBase), with_probability= prob_threshold, uncertainty_budget_mean=Γ, uncertainty_budget_variance=Γ)
    end

    # Chance constraints for the maximum output of generators:
    sum_deviations = sum([winderr[i] for i in 1:numfarms])
    for i in generatorlist
        addConstraint(m, pbar[i] - sum_deviations*alpha[i] >= buses[i].Pgmax, with_probability=gen_probability_threshold, uncertainty_budget_mean=Γ, uncertainty_budget_variance=Γ)
    end

    # Chance constraints on the upward (+) and downward(-) ramping:
    for (bus_idx, rampmax) in ramp
        addConstraint(m, -sum_deviations*alpha[bus_idx] >= rampmax, with_probability=gen_probability_threshold, uncertainty_budget_mean=Γ, uncertainty_budget_variance=Γ)
        addConstraint(m, -sum_deviations*alpha[bus_idx] <= -rampmax, with_probability=gen_probability_threshold, uncertainty_budget_mean=Γ, uncertainty_budget_variance=Γ)
    end


    # Set the objective function 
    @setObjective(m, Min, sum{buses[i].pi1*pbar[i]+buses[i].pi2*(pbar[i]^2 + sumvar*alpha[i]^2), i in generatorlist} + VoLL*sum(slack_bus))



    tic() # track the execution time
    status = solvechance(m,method=:Cuts, debug=false) # Solve the model
    solvetime = toq() # record the execution time
    if status != :Optimal
        # if the model isn't solved optimally, return NaNs
        return status, NaN, solvetime, fill(NaN, numbuses), fill(NaN, numbuses), fill(NaN, numlines), fill(NaN, numbuses-1)
    end
    # Return the optimal values if the model is solved optimally
    return status, getObjectiveValue(m), solvetime, getValue(alpha), getValue(pbar), getValue(barf), getValue(slack_bus)

end

# Read the input data:
logfilename, casefilename, windfilename, costsfilename, refbus, line_probability_threshold, gen_probability_threshold, linebufferamt, uniformalphas, mvaBase, loadscale, thermalLimitscale, extras = readconfig(ARGS[1])

# Read specific data for the test case:
numbuses, numgens, generatorlist, numbranches, buses, lines  = readcase(casefilename, costsfilename, refbus, loadscale, thermalLimitscale)

# Initialize index of specific (hydro, nucs, coal, gas) generators and transmission lines (loaded)
hydro_idx = int(readcsv(extras["hydro_indices"])[:])
nuclear_idx = int(readcsv(extras["nuclear_indices"])[:])
coal_idx = int(readcsv(extras["coal_indices"])[:])
gas_idx = int(readcsv(extras["gas_indices"])[:])
loaded_lines_idx = int(readcsv(extras["loaded_lines_file"])[:])
# Threshold for the loaded line:
loaded_lines_probability_threshold = float(extras["loaded_lines_probability_threshold"])

# Make sure that all generators are exactly one of the above types
@assert issubset(hydro_idx, generatorlist)
@assert issubset(nuclear_idx, generatorlist)
@assert issubset(coal_idx, generatorlist)
@assert issubset(gas_idx, generatorlist)

# Define the set of indices for hydro generators. Will be used to enforce the uniform participation factors.
fixed_alpha_subsets = [hydro_idx]

# Extract the names of additional input files from the outputs of redconfig.jl (line 142):
wind_forecast_file = extras["wind_forecast_file"]
wind_sigma_file = extras["wind_sigma_file"]
load_file = extras["load_file"]
hydro_limit_file = extras["hydro_limit_file"]
ramp_file = extras["ramp_file"]

# Read the input files:
wind_forecast = readcsv(wind_forecast_file)
wind_sigma = readcsv(wind_sigma_file)
loads = readcsv(load_file)
hydro_limit = readcsv(hydro_limit_file)
ramp_csv = readcsv(ramp_file)
ramp = [(int(ramp_csv[i,1]),ramp_csv[i,2]) for i in 1:size(ramp_csv,1)]

# Extract the name of the  out file from the outputs of redconfig.jl (line 142):
output_matfile = extras["output_matfile"]

# Specify the name of the files, from which the ranges on the mean and sigma should be read:
robust_mean_lower_file = extras["robust_mean_lower_file"]
robust_mean_upper_file = extras["robust_mean_upper_file"]
robust_sigma_lower_file = extras["robust_sigma_lower_file"]
robust_sigma_upper_file = extras["robust_sigma_upper_file"]
# Read the input range on the means and sigma 
robust_mean_lower = readcsv(robust_mean_lower_file)
robust_mean_upper = readcsv(robust_mean_upper_file)
robust_sigma_lower = readcsv(robust_sigma_lower_file)
robust_sigma_upper = readcsv(robust_sigma_upper_file)

# Set the value of the budget of uncertainty from the input file
robust_budget = extras["robust_budget"]

# determine if the robust ccopf or ccopf to be solved
ROBUST = false
try
    robust_budget = float(robust_budget)
    ROBUST = true
catch
    robust_budget = 0.0
end

# Display the name of the model to be solved
if ROBUST
    println("ROBUST CCOPF with budget $robust_budget")
else
    println("CCOPF")
end

# Time steps to be solved:
#numTimeSteps = size(wind_forecast,2)-1 # uncomment this line and comment the line below in order to run the full experiment. Note that this may take many hours to complete.
numTimeSteps = 10
println("$numTimeSteps time steps")

# Set up wind farms, using data from time step 1
@assert size(wind_forecast,1) == size(wind_sigma,1)
numfarms = size(wind_forecast,1)

farms = Farm[]
for i in 1:numfarms
    busidx = int(wind_forecast[i,1])
    f = Farm(busidx, 0, 0)
    f.colbarinv = Array(Float64, numbuses)
    push!(farms, f)
    setfarm(buses[busidx],i)
end

# Factorization, this only needs to be done once
B = genB(buses, lines)
Bhat = remove_col_and_row(B,refbus)
Bhatsp = sparse(Bhat)
Blu = lufact(Bhatsp)

# Compute terms of the inverse vector
rhs = zeros(numbuses-1)
for f in farms
    rhs[f.node] = 1
    f.colbarinv[1:numbuses .!= refbus] = Blu\rhs
    f.colbarinv[refbus] = 0
    rhs[f.node] = 0
end

# Arrays for storing solution values
objvals = Float64[]
solvetimes = Float64[]
statuses = ASCIIString[]
alphavals = {}
pbarvals = {}
barfvals = {}
slackvals = {}

for t in 1:numTimeSteps
    println("\n\n######## T = $t\n")
    for i in 1:size(loads,1)
        busidx = int(loads[i,1])
        Pd = float(loads[i,t+1])
        buses[busidx].Pd = Pd
    end
    for i in 1:numfarms
        busidx = int(wind_forecast[i,1])
        mean = loadscale*float(wind_forecast[i,t+1])
        
        std = loadscale*float(wind_sigma[i,t+1])
        farms[i].mean = mean
        farms[i].stddev = std
        if ROBUST
            #Scale the upper and lower bounds on the means and sigma ranges:
            mean_lb = loadscale*robust_mean_lower[i,t+1]
            mean_ub = loadscale*robust_mean_upper[i,t+1]
            std_lb = loadscale*robust_sigma_lower[i,t+1]
            std_ub = loadscale*robust_sigma_upper[i,t+1]
            if isnan(mean_lb) || isnan(mean_ub) || isnan(std_lb) || isnan(std_ub)
                error("Invalid data for farm at bus $busidx, T = $t")
            end
            # Write the mean and variance interval data for each wind farm:
            farms[i].deviation_mean_interval = (mean_lb-mean, mean_ub-mean)
            farms[i].deviation_sigmasq_interval = (std_lb^2,std_ub^2)
        else
            # if the ccopf to be solved, the mean is set to 0 and the variance is calculated accordingly:
            farms[i].deviation_mean_interval = 0.0
            farms[i].deviation_sigmasq_interval = std^2
        end

    end
    for i in 1:numbuses
        setMeanWindMinusLoad(buses[i], farms)
    end
    
    sumvar = 0.0
    for f in farms
        sumvar += f.stddev^2
    end
    println("-> sumvar ", sumvar, " stddev ", sqrt(sumvar))

    barb = (Float64[b.meanwindminusload for b in buses])[1:numbuses .!= refbus]
    Binvb = zeros(numbuses)
    Binvb[1:numbuses .!= refbus] = Blu\barb


    # Solve the ccopf model
    status, objval, solvetime, alphaval, pbarval, barfval, slackval =
      solve_ccopf(logfilename, refbus, line_probability_threshold, gen_probability_threshold, linebufferamt, mvaBase, loadscale, thermalLimitscale, buses, lines, farms, generatorlist, Binvb, Bhatsp, sumvar, hydro_idx, nuclear_idx, hydro_limit[:,t+1], ramp, robust_budget, loaded_lines_idx, loaded_lines_probability_threshold)

    # Write the optimal solution:  
    push!(objvals, objval)
    push!(solvetimes, solvetime)
    push!(statuses, string(status))
    push!(alphavals, alphaval)
    push!(pbarvals, pbarval)
    push!(barfvals, barfval)
    push!(slackvals, slackval)
end

# Prepare the output arrays (first column is generator index)
alphamat = Array(Float64, length(generatorlist), numTimeSteps+1)
pbarmat = Array(Float64, length(generatorlist), numTimeSteps+1)
barfmat = Array(Float64, length(lines), numTimeSteps+1)
slackmat = Array(Float64, length(buses)-1, numTimeSteps+1)

alphamat[:,1] = generatorlist
pbarmat[:,1] = generatorlist
barfmat[:,1] = 1:length(lines)
slackmat[:,1] = 1:length(buses)-1

# Create the output arrays
for t in 1:numTimeSteps
    for i in 1:length(generatorlist)
        alphamat[i,t+1] = alphavals[t][generatorlist[i]]
        pbarmat[i,t+1] = pbarvals[t][generatorlist[i]]
    end
    for i in 1:length(lines)
        barfmat[i,t+1] = barfvals[t][i]
    end
    for i in 1:(length(buses)-1)
        slackmat[i,t+1] = slackvals[t][i]
    end
end

#Write the output arrays into the output mat file
mfile = matopen(output_matfile, "w")
write(mfile, "objvals", objvals)
write(mfile, "solvetimes", solvetimes)
write(mfile, "status", statuses)
write(mfile, "alphavals", alphamat)
write(mfile, "pbarvals", pbarmat)
write(mfile, "barfvals", barfmat)
write(mfile, "slackvals", slackmat)
close(mfile)






