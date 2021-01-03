"""
    DispersalKernel{Radius}(formulation::F, kernel::K, cellsize::C, distancemethod::D)
    DispersalKernel{Radius}(; formulation=ExponentialKernel(1.0), cellsize=1.0, distancemethod=CentroidToCentroid())

Preferably use the keyword constructor to build the array from
a dispersal kernel function.
"""
struct DispersalKernel{R,N,K,F,C,D} <: AbstractKernel{R}
    "Neighborhood object"
    neighborhood::N
    "Kernel"
    kernel::K
    "Kernel formulation object"
    formulation::F
    "Simulation cell size"
    cellsize::C
    "Method for calculating distance between cells"
    distancemethod::D
end
DispersalKernel{R}(hood, kernel, formulation, cellsize, distancemethod) where R = begin
    # Build the kernel matrix
    newkernel = scale(buildkernel(formulation, distancemethod, cellsize, R))
    # Convert the kernel matrix to the type of the init array
    S = 2R + 1
    DispersalKernel{
        R,map(typeof, (hood, newkernel, formulation, cellsize, distancemethod))...
    }(
        hood, newkernel, formulation, cellsize, distancemethod
    )
end
DispersalKernel{R}(;
    hood=Window{R}(),            
    kernel=nothing,
    formulation=ExponentialKernel(),
    cellsize=1.0,
    distancemethod=CentroidToCentroid(),
) where R = begin
    DispersalKernel{R}(kernel, _buffer, formulation, cellsize, distancemethod)
end
ConstructionBase.constructorof(::Type{<:DispersalKernel{R}}) where R = DispersalKernel{R}

function DynamicGrids._setbuffer(n::DispersalKernel{R,<:Any,K,F,C,D}, _buffer) where {R,N,K,F,C,D} 
    hood = _setbuffer(neighborhood(n), buf)
    DispersalKernel{R,typeof(hood),K,F,C,D}(
        hood, n.kernel, n.formulation, n.cellsize, n.distancemethod
    )
end

DynamicGrids.radius(hood::DispersalKernel{R}) where R = R

DynamicGrids._buffer(hood::DispersalKernel) = hood._buffer
DynamicGrids.kernel(hood::DispersalKernel) = hood.kernel
cellsize(hood::DispersalKernel) = hood.cellsize
distancemethod(hood::DispersalKernel) = hood.distancemethod
formulation(hood::DispersalKernel) = hood.formulation

buildkernel(formulation, distancemethod, cellsize, r) = begin
    # The radius doesn't include the center cell, so add it
    sze = 2r + 1
    kernel = zeros(typeof(cellsize), sze, sze)
    r1 = r + one(r)
    # Paper: l. 97
    for x = 0:r, y = 0:r
        # Calculate the distance effect from the center cell to this cell
        prob = dispersalprob(formulation, distancemethod, x, y, cellsize)
        # Update the kernel value based on the formulation and distance
        kernel[x + r1, y + r1] = prob
        kernel[x + r1, -y + r1] = prob
        kernel[-x + r1, y + r1] = prob
        kernel[-x + r1, -y + r1] = prob
    end
    SMatrix{sze,sze}(kernel)
end

scale(x) = x ./ sum(x)


# build_dispersal_kernel(formulation, distancemethod, cellsize, r) = begin
#     # The radius doesn't include the center cell, so add it
#     sze = 2r + 1
#     kernels = [zeros(Float64, sze, sze) for thread in 1:Threads.nthreads()]
#     r1 = r + one(r)

#     # Calculate one quadrant
#     Threads.@threads for x = 0:r
#         kernel = kernels[Threads.threadid()]
#         for y = 0:r
#             kernel[x + r1, y + r1] = dispersalprob(formulation, distancemethod, x, y, cellsize)
#         end
#     end
#     kernel = sum(kernels)
#     println()
#     println(length(kernels))
#     display(kernels)
#     println()
#     # Update the other 3 quadrants
#     for x = 0:r, y = 0:r
#         # Just sum the threads, most are zero
#         prob = kernel[x + r1, y + r1]
#         kernel[x + r1, -y + r1] = prob
#         kernel[-x + r1, y + r1] = prob
#         kernel[-x + r1, -y + r1] = prob
#     end
#     # Normalise
#     kernel ./= sum(kernel)
# end


"""
Abstract supertype for methods of calculating distances and dispersal probabilities
between cells in a grid.
"""
abstract type DistanceMethod end

subsampling(method::DistanceMethod) = method.subsampling

"""
    CentroidToCentroid()

Calculates probability of dispersal between source and destination cell centroids
This is the obvious, naive method, but it will not handle low grid resolution well.
"""
struct CentroidToCentroid <: DistanceMethod end

dispersalprob(f, ::CentroidToCentroid, x, y, cellsize) =
    sqrt(x^2 + y^2) * cellsize |> f

# """
#     CentroidToCentroid(subsampling)
#     CentroidToCentroid(; subsampling=10.0)
#
# Calculates probability of dispersal between source cell centroid and destination cell area.
# """
# @columns struct CentroidToArea <: DistanceMethod
    # Field        | Default | Flat | Bounds      | Description
    # subsampling::Int | 10.0    | true | (2.0, 40.0) | "Subsampling for brute-force integration"
# end

# dispersalprob(f, dm::CentroidToArea, x, y, cellsize) = error("not implemented yet")

"""
    AreaToCentroid(subsampling)
    AreaToCentroid(; subsampling=10.0)

Calculates probability of dispersal between source cell area and destination centroid.
"""
Base.@kwdef struct AreaToCentroid{SS<:Number} <: DistanceMethod
    "Subsampling for numerical integration"
    subsampling::SS = Param(10.0; bounds=(2.0, 40.0))
end

@inline dispersalprob(f, dm::AreaToCentroid, x, y, cellsize) = begin
    prob = zero(cellsize)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, round(Int, subsampling(dm)))
    for j in range, i in range
        prob += sqrt((x + i)^2 + (y + j)^2) * cellsize |> f
    end
    prob / subsampling(dm)^2
end

"""
    AreaToArea(subsampling)
    AreaToArea(; subsampling=10.0)


Calculates probability of dispersal between source and destination cell areas.
"""
Base.@kwdef struct AreaToArea{SS<:Number} <: DistanceMethod
    subsampling::SS = Param(10.0; bounds=(2.0, 40.0))
end

@inline dispersalprob(f, dm::AreaToArea, x, y, cellsize) = begin
    prob = zero(cellsize)
    # Get the center point of the first cell (for both dimensions)
    centerfirst = 1 / subsampling(dm) / 2 - 0.5
    centerlast = centerfirst * -1
    range = LinRange(centerfirst, centerlast, round(Int, subsampling(dm)))
    for i in range, j in range
        for a in range, b in range
            prob += sqrt((x + i + a)^2 + (y + j + b)^2) * cellsize |> f
        end
    end
    prob / subsampling(dm)^4
end


"""
Abstract supertype for functors that calculate the probability of
dispersal between two points.

Concrete implementations must define functor methods with the form:

```julia
(k::SomeKernel)(distance) = # do something with `distance` and `k`"
```

Using an anonymous funciton for this would not give you rebuildable
model parameters.
"""
abstract type KernelFormulation end

"""
    ExponentialKernel(λ)

Probability of dispersal with a negatitve exponential relationship to distance.
"""
@Base.kwdef struct ExponentialKernel{P} <: KernelFormulation
    "Parameter for adjusting spread of dispersal propability"
    λ::P = Param(1.0, bounds=(0.0, 2.0))
end
(f::ExponentialKernel)(distance) = exp(-distance / f.λ)
