using CellularAutomataBase: normalizeframe

" Parametrizer to use with Optim.jl or similar "
struct Parametriser{R,OB,L,NR,TS,OP}
    ruleset::R
    objective::OB
    loss::L
    nreplicates::NR
    tstop::TS
    output::OP
end
Parametriser(ruleset, objective, loss, nreplicates, tstop) = begin
    output = ArrayOutput(ruleset.init, tstop)
    Parametriser(ruleset, objective, loss, nreplicates, tstop, output)
end

" Objective function for the parametriser "
(p::Parametriser)(params) = begin
    # Rebuild the rules with the current parameters
    names = fieldnameflatten(p.ruleset.rules, Real)
    println("Parameters: ", collect(zip(names, params)))
    p.ruleset.rules = Flatten.reconstruct(p.ruleset.rules, params, Real)
    cumsum = @distributed (+) for i = 1:p.nreplicates
        output = deepcopy(p.output)
        sim!(output, p.ruleset; tstop = p.tstop)
        prediction = outputtoprediction(p.objective, output)
        loss = value(p.loss, target(p.objective), prediction, AggMode.Sum())
        println("replicate: ", i, " - loss: ", loss)
        loss
    end
    meanloss = cumsum ./ p.nreplicates
    println("mean loss: ", meanloss, "\n")
    return meanloss
end


"""
AbstractObjectives map simulation outputs to predictions that 
can be compared to target data using a loss function.

THey must implement `outputtoprediction` and `targets` methods.
"""
abstract type AbstractObjective end


"""
    outputtoprediction(obj::AbstractObjective, output::AbstractOutput)
Methods that map an objective object and a simulation output to a 
prediction array.
"""
function outputtoprediction end

"""
    target(obj::AbstractObjective)
Returns a target array given an AbstractObjective. The target must match the size and 
dimensions of the prediction array returned by `outputtoprediction`.
"""
function target end



"""
Implementation of a loss objective that converts cell data to regional
presence/absence and compares to a target of regional occurance data.
"""
struct RegionObjective{DT,RL,OC,FS} <: AbstractObjective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
end

target(obj::RegionObjective) = obj.occurance

outputtoprediction(obj::RegionObjective, output) = begin
    regions, steps = size(obj.occurance)
    frames = length(output)
    # Allocate arrays for steps and set all cells to zero 
    outputsteps = [similar(output[1]) for f in 1:steps]
    fill!.(outputsteps, zero(eltype(outputsteps[1])))
    # Get the mean population for steps from the frames in each step
    for frame in 1:frames
        step = stepfromframe(obj.framesperstep, frame)
        outputsteps[step] .+= output[frame]
    end
    # Divide all cells by frames per step to get the mean population
    map(s -> s ./= obj.framesperstep, outputsteps)

    # Allocate a boolean array to contain our presence/absence predictions
    prediction = zeros(Bool, size(obj.occurance))
    # Convert mean cell populations to regional prescence/absence
    for t in 1:steps
        for r in 1:regions
            prediction[r, t] = (sum((obj.regionlookup .== r) .& (outputsteps[t] .> 0)) ./
                       sum((obj.regionlookup .== r))) > obj.detectionthreshold
        end
    end
    prediction 
end

stepfromframe(framesperstep, t) = (t - one(t)) ÷ framesperstep + one(t)


"""
An image procesor to visualise the model fit, for a live version of
the region fitting optimiser.

Fields:
`objective` : a RegionObjective object
`truepositivecolor` : color of true positive fit, etc.
`falsepositivecolor`
`truenegativecolor`
`falsenegativecolor`
`maskcolor` : color when a cell region of zero or lower
"""
struct ColorRegionFit{O<:RegionObjective,TP,FP,TN,FN,M} <: AbstractFrameProcessor
    objective::O
    truepositivecolor::TP
    falsepositivecolor::FP
    truenegativecolor::TN
    falsenegativecolor::FN
    maskcolor::M
end

CellularAutomataBase.frametoimage(p::ColorRegionFit, output, frame, t) = begin
    step = stepfromframe(p.objective.framesperstep, t)
    img = similar(frame, RGB24)
    for i in CartesianIndices(frame)
        region = p.objective.regionlookup[i]
        img[i] = if region > zero(region)
            x = frame[i]
            if p.objective.occurance[region, step]
                x == zero(x) ? rgb(p.falsenegativecolor) : rgb((x .* p.truepositivecolor))
            else
                x == zero(x) ? rgb(p.truenegativecolor) : rgb((x .* p.falsepositivecolor))
            end
        else
           rgb(p.maskcolor)
        end
    end
    img
end

rgb(c::RGB24) = c
rgb(c::Tuple) = RGB24(c...)
rgb(c::Number) = RGB24(c)
