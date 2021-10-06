module RootsByDistribution

using DataStructures, Roots, Sobol
import Base.Iterators

export roots, bracket_roots

# Find n roots of the continuous function f defined on the image of [0,1] under
# the given transform.
#
# First, brackets are found using bracket_roots. Then, Roots.find_zero is used
# to find a root in each bracket.
#
# The transform determines the distribution from which points are sampled,
# which can influence the number of necessary evaluations, see the
# bracket_roots.
#
# This function will not terminate if f has less than n brackets in the image
# of [0,1] under the given transform.

function roots(f, n, transform=Float64)
    map(bracket_roots(f, n, transform)) do (lo, hi)
        lo == hi ? lo : find_zero(f, (lo, hi))
    end
end


# Find n intervals containing roots of the continuous function f by sampling
# points from a low-discrepancy sequence in [0,1] mapped under the given
# transform. The sequence contains the endpoints 0 and 1.
#
# More precisely, the sampled points are of the form transform(x_i) where x₀=0,
# x₁=1 and xᵢ for i≥2 is the Sobol sequence in [0,1].
#
# This function will not terminate if f has less than n brackets in the image
# of [0,1] under transform.

function bracket_roots(f, n, transform=Float64)
    nt = (typeof ∘ transform ∘ Float64)(0)
    intervals = Vector{Tuple{nt, nt}}()
    n <= 0 && return intervals

    iter(seq) = Iterators.map(seq) do z
        x = transform(z)
        s = (sign ∘ f)(x)
        (x, s)
    end

    # The Sobol sequence does not sample the endpoints, but we want to see
    # whether they are roots.

    seq0 = (0.0, 1.0)
    (x₀, s₀), (x₁, s₁) = iter(seq0)
    d = SortedDict((x₀ => s₀, x₁ => s₁))
    # k is the number of sign changes detected.
    k = count((s₀ == 0,
               s₁ == 0,
               s₀ == -s₁ && s₀ != 0))

    seq1 = Iterators.map(first, SobolSeq(1))
    for (x, s) in iter(seq1)
        k >= n && break
        xlast, slast = deref((d, searchsortedlast(d, x)))
        xfirst, sfirst = deref((d, searchsortedfirst(d, x)))
        (x == xlast || x == xfirst) && continue
        push!(d, x => s)
        slast == -sfirst && slast != 0 && (k = k - 1)
        k = k + count((s == 0,
                       s != 0 && s == -slast,
                       s != 0 && s == -sfirst))
    end

    skip = false
    for ((x, s), (y, t)) in zip(d, Iterators.drop(d, 1))
        !skip && s == 0 && push!(intervals, (x, x))
        skip = false
        t == 0 && (push!(intervals, (y, y)); skip = true)
        s != 0 && s == -t && push!(intervals, (x, y))
    end

    intervals
end

end
