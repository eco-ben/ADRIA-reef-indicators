using Distances
using Clustering
using FLoops
using ADRIA
using JuliennedArrays

"""
    _complexity(x::AbstractMatrix{<:Real})::AbstractMatrix{Float64}

Compute Complexity (CE) of an Matrix `x` of shape \$T ⋅ S\$, where \$T\$ is total number of
time steps and \$S\$ is number of scenarios.

# Arguments
- `x` : series matrix of shape \$T ⋅ S\$

# Returns
Vector of \$N\$ elements
"""
function _complexity(x::AbstractMatrix{<:Real})::Vector{Float64}
    return vec(sqrt.(sum(diff(Matrix(x); dims=1) .^ 2; dims=1)) .+ 1)
end

"""
    correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}

Compute Correction Factor (CF) between two time series complexities `ce_i` and `ce_j`.

# Arguments
- `ce_i` : Time series `i`
- `ce_j` : Time series `j`

# Returns
Float64

# Examples
```julia
julia> ce = complexity([[1, 2, 3] [1, 3, 4]])
julia> correction_factor(ce[1], ce[2])
Float64:
 2.5
 ```
"""
function correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}
    return max(ce_i, ce_j) / min(ce_i, ce_j)
end

function complexity_invariance_distance(
    data::AbstractMatrix{<:Real};
    distance=:euclidean
)::AbstractMatrix{Float64}
    # Compute complexity vector
    complexity = _complexity(data)

    # Create empty Matrix
    data_size = size(data, 2)
    cid_matrix::AbstractMatrix{Float64} = zeros(data_size, data_size)
    dist_complexity = [(data[:, i], complexity[i]) for i in 1:data_size]

    local weights::Vector{Float64}
    if distance == :weuclidean
        # [1, 1/2, 1/3, ..., 1/n]
        weights = sqrt.(1 ./ (1:size(data, 1)))
        #cid_matrix = create_distance_matrix(dist_complexity; dist_fn = weuclidean(x, y, weights))
        
        return cid_matrix
    end
    dist_fn(x, y) = (distance==:euclidean) ? euclidean(x, y) : weuclidean(x, y, weights)
    
    for ii in axes(cid_matrix, 1)
        Threads.@threads for jj in axes(cid_matrix, 2)
            if ii == jj || !iszero(cid_matrix[ii, jj])
                continue
            end

            @views cid_matrix[ii, jj] = complexity_invariance(dist_complexity[ii], dist_complexity[jj], dist_fn)
        end

        cid_matrix[:, ii] .= cid_matrix[ii, :]
    end
    
    # cid_matrix = create_distance_matrix(dist_complexity)
    return cid_matrix
end

function complexity_invariance((data_x, complexity_x), (data_y, complexity_y), dist_fn)
    ed = dist_fn(data_x, data_y)
    cf = correction_factor(complexity_x, complexity_y)
    return ed * cf
end

"""
    cluster_series(data::AbstractMatrix{<:Real}, n_clusters::Int64, method::Symbol=:kmedoids, distance::Symbol=:euclidean)::Vector{Int64}

Hierarchically cluster \$S\$ scenarios with \$T\$ time steps each.

# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is
  number of scenarios
- `n_clusters` : Number of clusters determined _a priori_
- `method` : Clustering method. Defaults to `:kmedoids`
- `distance` : Switch between Euclidean (`:euclidean`) or weighted Euclidean (`:weuclidean`)
distance measurements. Defaults to `:euclidean`

# Returns
- Cluster ids indicating each scenario cluster assignment.

# References
1. Steinmann, P., Auping, W.L., Kwakkel, J.H., 2020.
   Behavior-based scenario discovery using time series clustering.
   Technological Forecasting and Social Change 156, 120052.
   https://doi.org/10.1016/j.techfore.2020.120052

2. Batista, G.E.A.P.A., Keogh, E.J., Tataw, O.M., de Souza, V.M.A., 2014.
   CID: an efficient complexity-invariant distance for time series.
   Data Min Knowl Disc 28, 634-669.
   https://doi.org/10.1007/s10618-013-0312-3
"""
function cluster_series(
    data::AbstractMatrix{<:Real},
    n_clusters::Int64;
    method::Symbol=:kmedoids,
    distance::Symbol=:euclidean
)::Vector{Int64}
    # Calculate distantes matrix
    distances = complexity_invariance_distance(data; distance=distance)

    if method == :kmedoids
        return kmedoids(distances, n_clusters).assignments
    end

    # Return hierarchical clustering with n_clusters
    dendogram = hclust(distances; linkage=:average)
    return cutree(dendogram; k=n_clusters)
end


"""
    distance_matrix!(reef_geoms::Vector, dist::Matrix{Float64}, func::Function)::Nothing

Create a distance matrix using the provided function.

# Arguments
- `reef_geoms` : Geometries to calculate distances between
- `dist` : distanec/adjacency matrix
- `func` : Distance calculation function
"""
function distance_matrix!(reef_geoms::Vector, dist::Matrix{Float64}, func::Function)::Nothing
    for ii in axes(dist, 1)
        Threads.@threads for jj in axes(dist, 2)
            if ii == jj || !iszero(dist[ii, jj])
                continue
            end

            @views dist[ii, jj] = func(reef_geoms[ii], reef_geoms[jj])
        end

        dist[:, ii] .= dist[ii, :]
    end

    return nothing
end

"""
    create_distance_matrix(reef_geoms::Vector)::Matrix

Calculate matrix of unique distances between reefs.

# Returns
Distance between reefs in meters
"""
function create_distance_matrix(reef_geoms::Vector)::Matrix{Float64}
    n_reefs = size(reef_geoms, 1)
    dist = zeros(n_reefs, n_reefs)
    distance_matrix!(reef_geoms, dist)

    return dist
end

# start with this one 
@inline function distance_matrix!(reef_geoms::Vector{Tuple{Vector{Float64}, Float64}}, dist::Matrix{Float64})
    return distance_matrix!(reef_geoms, dist, (x, y) -> complexity_invariance(x, y, dist_fn))
end

# @inline function distance_matrix!(reef_geoms::Vector{Tuple{Vector{Float64}, Float64}}, dist::Matrix{Float64}, dist_fn)
#     return distance_matrix!(reef_geoms, dist, complexity_invariance)
# end
