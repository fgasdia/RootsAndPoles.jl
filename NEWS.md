# RootsAndPoles.jl Release Notes

Important notices about MAJOR and MINOR version releases.

- `v1.5.0`: Minimum Julia version now `v1.7`. Removes custom implementation of fast `delaunayedges`.

- `v1.4.0`: Provides a new function for preparing data for plotting, `getplotdata`.
A fix for `rectangulardomain` is also included that, as a side effect, shifts the nodes produced by the old `rectangulardomain` by a half spacing.
This may cause slightly different root and pole locations from `grpf`.

- `v1.3.0`: `GRPFParams` no longer parameterizes the type of the `tolerance` field and instead it is fixed as `Float64`.
There is no reason `tolerance` would have been anything other than `Float64` because the precision is limited to `Float64` by the underlying `VoronoiDelaunay.jl` triangulation procedure.

- `v1.2.0`: `GRPFParams` are now `mutable` structs.
Although this can be considered a breaking change because two structs which are equal are now not identical, in the sense that `a == b` but `a !== b`, most users should not have any code broken by this change.

- `v1.1.0`: Calls to the user-provided complex function can be multithreaded using Julia's `@threads` capability.
By default, `multithreading = false`.
It can be enabled in the `GRPFParams` struct.
