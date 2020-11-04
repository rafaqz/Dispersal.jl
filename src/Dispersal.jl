module Dispersal
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Dispersal

using ConstructionBase,
      Colors,
      Dates,
      DimensionalData,
      Distributed,
      Distributions,
      DocStringExtensions,
      FieldDefaults,
      FieldDocTables,
      FieldMetadata,
      Flatten,
      LinearAlgebra,
      LossFunctions,
      Mixers,
      PoissonRandom,
      Reexport,
      Setfield,
      Statistics

@reexport using DynamicGrids

using LossFunctions: ZeroOneLoss
using DimensionalData: Time

import Base: getindex, setindex!, lastindex, size, length, push!


import DynamicGrids: applyrule, applyrule!, applyrule, applyrule!,
       neighbors, sumneighbors, neighborhood, setneighbor!, mapsetneighbor!,
       radius, gridsize, mask, overflow, cellsize, ruleset, isinbounds, inbounds, ismasked,
       extent, currenttime, currenttimestep, timestep, tspan,
       buffer, WritableGridData, SimData, aux, unwrap

import ConstructionBase: constructorof

import FieldMetadata: @default, @description, @bounds, @flattenable,
                      default, description, bounds, flattenable


export AbstractDispersalKernel, DispersalKernel

export KernelFormulation, ExponentialKernel

export DistanceMethod, CentroidToCentroid, CentroidToArea, AreaToArea, AreaToCentroid

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       SwitchedInwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal


export AlleeExtinction, JumpDispersal, HumanDispersal

export BatchGroups, HierarchicalGroups

export ExactExponentialGrowth, ExactLogisticGrowth

export MaskGrowthMap, ExactExponentialGrowthMap, ExactLogisticGrowthMap

export LayerCopy

export Parametriser, AbstractObjective, SimpleObjective, RegionObjective, RegionOutput,
       ColorRegionFit, Accuracy

export ThreadedReplicates, DistributedReplicates, SingleCoreReplicates

export targets, predictions


const FIELDDOCTABLE = FieldDocTable((:Description, :Default, :Bounds),
                                    (description, default, bounds);
                                    truncation=(100,40,100))

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

include("rules/mixins.jl")
include("downsampling.jl")
include("layers.jl")
include("rules/kernel/common.jl")
include("rules/kernel/inwards.jl")
include("rules/kernel/outwards.jl")
include("rules/growth.jl")
include("rules/human.jl")
include("rules/jump.jl")
include("rules/allee.jl")
include("optimisation/optimisation.jl")
include("optimisation/objectives.jl")
include("optimisation/output.jl")


end # module
