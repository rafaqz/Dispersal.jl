using DynamicGrids: normalise, minval, maxval, ismasked, rgb24

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
struct ColorRegionFit{O,P,N,TZ,FZ,M} <: GridProcessor
    objective::O
    truescheme::P
    falsescheme::N
    truezerocolor::TZ
    falsezerocolor::FZ
    maskcolor::M
end

"""
    grid2image(p::ColorRegionFit, output, ruleset, frame, t)

Visualise the match between predictions and observed regional occupancy 
during live simulations.
"""
DynamicGrids.grid2image(p::ColorRegionFit, output::ImageOutput, ruleset::Ruleset, frame, t) = begin
    step = stepfromframe(p.objective, t)
    img = similar(frame, RGB24)
    obj = p.objective
    min, max = minval(ruleset), maxval(ruleset)
    for i in CartesianIndices(frame)
        region = p.objective.regionlookup[i]
        img[i] = if !(p.maskcolor isa Nothing) && ismasked(mask(ruleset), i) 
            p.maskcolor
        elseif region > zero(region)
            x = frame[i]
            normed = normalise(x, min, max)
            if p.objective.occurance[region, step]
                if !(p.truezerocolor isa Nothing) && normed == zero(normed) 
                    rgb24(p.falsezerocolor)
                else
                    rgb24(p.truescheme, normed)
                end
            else
                if !(p.falsezerocolor isa Nothing) && normed == zero(normed) 
                    rgb24(p.truezerocolor)
                elseif x > obj.detectionthreshold 
                    rgb24(p.falsescheme, normed)
                else
                    rgb24(p.truescheme, normed)
                end
            end
        else
            p.maskcolor
        end
    end
    img
end