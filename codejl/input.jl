# Code for processing input data for CCOPF simulation
# We define Julia types to store the bus and line data.

using MAT

type Bus
    nodeID::Int
    kind::Int
    Pd::Float64
    Qd::Float64
    Pg::Float64
    Qg::Float64
    Pgmax::Float64
    Qgmax::Float64
    Pgmin::Float64
    Qgmin::Float64
    pi2::Float64
    pi1::Float64
    Farmids::Vector{Int}
    meanwindminusload::Float64
    qobjcoeff::Float64
    genids::Vector{Int}
    outlist::Vector{Int}
    inlist::Vector{Int}
    function Bus(nodeID, kind, Pd, Qd)
        b = new(nodeID, kind, Pd, Qd)
        b.Pg = 0
        b.Qg = 0
        b.Pgmax = 0
        b.Qgmax = 0
        b.pi1 = 0
        b.pi2 = 0
        b.Farmids = Int[]
        b.meanwindminusload = 0
        b.genids = Int[]
        b.outlist = Int[]
        b.inlist = Int[]
        return b
    end
end

# copied from bus.py
function setg(b::Bus, genidx, Pg, Qg, Pgmax, Pgmin)
    b.Pg += Pg
    b.Qg += Qg
    b.Pgmax += Pgmax
    b.Pgmin += Pgmin
    if b.kind == 1
        warn("Generator $genidx was assigned to bus $(b.nodeID), but this bus has type 1")
    end
    b.kind = 2
    push!(b.genids,genidx)
end

getLoad(b::Bus) = b.Pd
getPg(b::Bus) = b.Pg
isgen(b::Bus) = b.kind == 2
getb(b::Bus) = b.Pg - b.Pd
getLoadMinusMeanWind(b::Bus) = -b.meanwindminusload
function setfarm(b::Bus, farmid)
    println("node ", b.nodeID, " is farm ", farmid)
    push!(b.Farmids, farmid)
end
function setMeanWindMinusLoad(b::Bus, farms)
    b.meanwindminusload = -b.Pd
    for farmid in b.Farmids
        b.meanwindminusload += farms[farmid].mean
    end
end
setqobjcoeff(b::Bus, coeff) = (b.qobjcoeff = coeff)

type Line
    arcID::Int
    tail::Int # the "to" node
    head::Int # the "from" node
    y::Float64 # the susceptance value
    x::Float64 # the reactance value
    u::Float64 # the capacity of the line
    distance_scale::Float64 # this will be used to scale u
    Line(arcID, tail, head, y, x, u, d) = new(arcID, tail, head, y, x, u, d)
end

Line(arcID, tail, head, y, distance_scale) = Line(arcId, tail, head, y, 1/y, 0.0, distance_scale)

getThermalCapacity(l::Line, mvaBase) = l.u
getSyncCapacity(l::Line, mvaBase) = l.y*mvaBase

type Farm
    node::Int
    mean::Float64
    stddev::Float64
    colbarinv::Vector{Float64}
    deviation_mean_interval
    deviation_sigmasq_interval
    Farm(idx, mean, stddev) = new(idx, mean, stddev, Float64[], mean, stddev^2)
end


function readcase(casefilename, costsfilename, refbus, loadscale, thermalLimitScale)
    # takes a MATPOWER case file, written out to .mat format
    file = matopen(casefilename)
    case = read(file, "mpc")
    close(file)

    ## bus data
    #	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
    busmat = case["bus"]
    buses = Bus[]
    for i in 1:size(busmat,1)
        @assert busmat[i,1] == i # indexed in order
        bustype = int(busmat[i,2])
        Pd = busmat[i,3]
        Qd = busmat[i,4]
        Gs = busmat[i,5]
        Bs = busmat[i,6]
        area = busmat[i,7]
        Vm = busmat[i,8]
        Va = busmat[i,9]
        baseKV = busmat[i,10]
        zone = busmat[i,11]
        Vmax = busmat[i,12]
        Vmin = busmat[i,13]

        b = Bus(i, bustype, loadscale*Pd, Qd)
        push!(buses, b)
    end

    ## generator data
    #%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    generatorlist = Int[]
    genmat = case["gen"]
    for i in 1:size(genmat,1)
        busidx = int(genmat[i,1])
        Pg = genmat[i,2]
        Qg = genmat[i,3]
        Pgmax = loadscale*genmat[i,9]
        Pgmin = loadscale*genmat[i,10]
        push!(generatorlist, busidx)
        setg(buses[busidx],i,Pg, Qg, Pgmax, Pgmin)
    end

    for i in 1:length(buses)
        if buses[i].kind == 2 && length(buses[i].genids) == 0
            warn("Bus $i is listed as a generator, but no corresponding generator information found")
        end
    end


    ## branch data
    #	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
    branchmat = case["branch"]
    lines = Line[]
    for i in 1:size(branchmat,1)
        fbus = int(branchmat[i,1])
        tbus = int(branchmat[i,2])
        x = branchmat[i,4]
        y = 1/x
        u = branchmat[i,6]


        push!(buses[fbus].outlist, i)
        push!(buses[tbus].inlist, i)

        l = Line(i, tbus, fbus, y, x, u*thermalLimitScale, thermalLimitScale)
        push!(lines,l)
    end

    generatorlist = unique(generatorlist)


    if costsfilename != "none"
        costsmat = readcsv(costsfilename, Float64)
        for i in 1:size(costsmat,1)
            busid = int(costsmat[i,1])
            @assert length(buses[busid].genids) == 1
            buses[busid].pi2 = costsmat[i,2] # quadratic coefficient
            buses[busid].pi1 = costsmat[i,3] # linear coefficient
        end
    else
        for g in generatorlist
            buses[g].pi2 = 1 # dummy cost
            buses[g].pi1 = 0
        end
    end

    return length(buses), 1:size(genmat,1), generatorlist, length(lines), buses, lines
end


function readconfig(configfilename)
    println("\nreading config $configfilename")
    refbus = 0
    uniformalphas = false
    
    lines = readlines(open(configfilename,"r"))

    numlines = length(lines)

    windfilename = "NONE"
    costsfilename = "none"
    logfilename = "NONE"
    casefilename = "NONE"

    lookingforend = true
    line_probability_threshold = 0.005
    mvaBase = 100
    gen_probability_threshold = 0.05
    linebufferamt = 0
    loadscale = thermalLimitscale = 1.0

    extras = Dict()

    for l in lines
        beginswith(l,'#') && continue
        
        thisline = split(l)
        length(thisline) > 0 || continue
        if thisline[1] == "END"
            break
        elseif thisline[1] == "case"
            casefilename = thisline[2]
        elseif thisline[1] == "wind"
            windfilename = thisline[2]
        elseif thisline[1] == "costs"
            costsfilename = thisline[2]
        elseif thisline[1] == "refbus"
            refbus = int(thisline[2])
        elseif thisline[1] == "line_probability_threshold"
            line_probability_threshold = float(thisline[2])
            println(">>>> line_probability_threshold = $line_probability_threshold")
        elseif thisline[1] == "gen_probability_threshold"
            gen_probability_threshold = float(thisline[2])
            println(">>>> gen_probability_threshold = $gen_probability_threshold")
        elseif thisline[1] == "linebufferamt"
            linebufferamt = float(thisline[2])
            println(">>>> linebufferamt = $linebufferamt")
        elseif thisline[1] == "uniformalphas"
            uniformalphas = true
            println(">>>> uniform alphas")
        elseif thisline[1] == "mvaBase"
            mvaBase = float(thisline[2])
            println(">>>> mvaBase = $mvaBase")
        elseif thisline[1] == "loadscale"
            loadscale = float(thisline[2])
            println(">>>> loadscale = $loadscale")
        elseif thisline[1] == "thermalLimitscale"
            thermalLimitscale = float(thisline[2])
            println(">>>> thermalLimitscale = $thermalLimitscale")
        elseif thisline[1] == "logfile"
            logfilename = thisline[2]
            println(">>>> logfilename = $logfilename")
        elseif thisline[1] == "factor" || thisline[1] == "rmatrix" || thisline[1] == "line_sync_nu"
            println(">>>> ignoring $(thisline[1])")
        else
            extras[thisline[1]] = thisline[2]
            println(">>>> $(thisline[1]) = $(thisline[2])")
        end
    end

    return logfilename, casefilename, windfilename, costsfilename, refbus, line_probability_threshold, gen_probability_threshold, linebufferamt, uniformalphas, mvaBase, loadscale, thermalLimitscale, extras
end


# includes farmstobuses and setMeanWindMinusLoad
function readwind(windfilename, buses, loadscale)
    println("reading wind file $windfilename")

    lines = readlines(open(windfilename,"r"))

    numfarms = int(split(lines[1])[1])

    farms = Farm[]

    for i in 1:numfarms
        id, mean, std = split(lines[i+1])
        push!(farms, Farm(int(id), loadscale*float(mean), loadscale*float(std)))
        println("farm $i : $id $mean $std")
    end
    
    for i in 1:numfarms
        setfarm(buses[farms[i].node],i)
    end
    for b in buses
        setMeanWindMinusLoad(b, farms)
    end

    return numfarms, farms
end
