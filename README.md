# RootsByDistribution.jl

A Julia package for finding roots of a continuous function when the approximate
distribution of the roots is known.

The algorithm finds brackets (root-isolating intervals) by sampling points from
a low-discrepancy sequence adapted to the given distribution. The [Sobol
package](https://github.com/stevengj/Sobol.jl) is used for this purpose.

The roots are then found using a bisection-like algorithm from the [Roots
package](https://github.com/JuliaMath/Roots.jl).
