# Gather/scatter analysis (2.5D_inplane)

This is the analysis of the field gather (`FIELD_get_elem_sub`) and scatter
(`FIELD_add_elem`) operations that the earlier profiling flagged as about
34% of runtime, the next bottleneck after the force apply. All numbers are
on the 2.5D_inplane case (nelem 12800, ngll 5, npoin 205761, ndof 2,
elastic isotropic).

## What the two operations actually do

Global fields (`displ`, `veloc`, `accel`) are stored as `(npoin, ndof)`.
Fortran is column major, so for a node k the two components `F(k,1)` and
`F(k,2)` are `npoin` apart in memory, about 1.6 MB apart here. Each field
is about 3.3 MB.

Gather (per element, per node):

    k = ibool(i,j,e)
    fout(i,j,:) = Fin(k,:)      ! reads Fin(k,1), Fin(k,2), 1.6 MB apart

Scatter (per element, per node):

    k = ibool(i,j,e)
    Fout(k,:) = Fout(k,:) + fin(i,j,:)   ! read-modify-write, same stride

The indirection through `ibool` is inherent to the assembly, so the
question is how much the memory layout and access pattern cost, and what
can be changed.

## Method

I read the real `ibool` from `ibool_sem2d.dat` and microbenchmarked the two
inner loops directly, comparing the current `(npoin,ndof)` layout against a
`(ndof,npoin)` layout where a node's two components sit next to each other
(one cache line instead of two). I ran each in RCM element order (what the
solver actually uses, `OPT_RENUMBER=.true.`) and in a random element order
(to probe the poor-locality, memory-bound regime). Per-node timings, median
of three runs:

| Regime | gather (npoin,ndof -> ndof,npoin) | scatter (npoin,ndof -> ndof,npoin) |
|---|---|---|
| RCM order (what the solver does) | 0.68 -> 0.68 ns (1.00x) | 0.90 -> 0.80 ns (1.13x) |
| Random order (poor locality) | 1.30 -> 1.23 ns (1.05x) | 1.87 -> 1.28 ns (1.45x) |

## Finding 1: the data-layout change is not the win it looks like on paper

The cache-line argument says `(ndof,npoin)` should roughly halve the memory
traffic per node, so up to 2x. It does not play out. In the order the
solver actually runs (RCM), the layout change does nothing for gather and
about 1.13x for scatter. Even in the pessimistic random-order regime it is
1.05x gather and 1.45x scatter.

The reason gather barely moves is that at 0.68 ns per node it is not waiting
on memory, it is bound by the instruction stream (the indirect load of
`ibool`, the address arithmetic, the writes into `fout`). The two field
components being on separate cache lines does not matter when both lines are
already resident. Scatter benefits a bit more because it is a
read-modify-write, so halving the lines touched saves real traffic,
especially when locality is poor.

Weighting by how often each runs (two gathers of d and v, one scatter per
element), the `(ndof,npoin)` transpose works out to roughly 1.5 to 2% end
to end in the realistic regime. That is a small payoff for an invasive
change: the layout of `displ`/`veloc`/`accel` is assumed `(npoin,ndof)`
across time integration, boundary conditions, sources, receivers, energy,
and output, so transposing it touches a lot of code with real correctness
risk. Based on this measurement I do not think the transpose is worth it as
the first move.

## Finding 2: the RCM element ordering is already doing real work

Going from RCM order to random order roughly doubles both gather (0.68 ->
1.30 ns) and scatter (0.90 -> 1.87 ns). So the existing reverse
Cuthill-McKee element renumbering is already buying about a 2x locality
improvement for free. That also means there is not much more to extract
from node reordering on top of it, the easy locality is already captured.

## Finding 3: the real win is not gathering `veloc` when it is not used

`compute_Fint` gathers both `d` and `v` for every element every timestep
and passes both to `MAT_Fint`. But `MAT_Fint` only reads `v` on the
Kelvin-Voigt path (`MAT_KV_add_etav`). For a plain elastic material, which
is this whole test case and a large fraction of real runs, `v` is gathered
and never touched.

I tested skipping that gather (a one line guard) and measured with the
`OPT_FINT_PROFILE` wall timer:

| Build | compute_Fint wall | time loop wall |
|---|---|---|
| gathers v (current) | 18.81 s | 27.57 s |
| skips unused v | 16.42 s | 24.54 s |

That is a 12.7% drop in `compute_Fint`, which is about 9% of total runtime,
and the seismograms come out bit for bit identical (0 of 1610 values
differ), as they must, since `v` was dead for this material. This is a
bigger effect than the layout transpose and it is nearly free and low risk.

The catch is that it has to be conditional. Kelvin-Voigt uses `v`, so the
gather cannot just be deleted, it has to be gated on whether the material
actually reads `v`.

## Implemented

I implemented the conditional gather. To keep `solver.f90` from hard coding
which materials read `v`, I added a small query `MAT_needs_veloc(matpro)` in
`mat_gen`, next to `MAT_Fint` where that knowledge belongs. It returns
`MAT_isKelvinVoigt(matpro)` today, and the comment says to extend it if any
other material path starts reading `v`. `compute_Fint` (both the serial and
the OpenMP path) now gathers `v` only when `MAT_needs_veloc` is true.

Verified:

- serial output is bit for bit identical to the original baseline (0 of
  1610 values differ in Ux and Uz), as expected since `v` was dead for this
  elastic material.
- the OpenMP path at 4 threads gives the same round-off level agreement it
  already had from the coloring (max abs diff 1.14e-13 on Ux, 1.39e-17 on
  Uz), with no new differences introduced by the gather change.
- `compute_Fint` wall time drops from 18.81 s to 16.42 s, the 12.7% (about
  9% end to end) measured above.

## Recommendation for what is next

1. (Done) The conditional `v` gather, about 9% end to end here, bit
   identical results, small and local change, no data structure churn. It
   also stacks with the OpenMP work, since it shrinks the parallel region's
   memory traffic.
2. Hold off on the `(npoin,ndof)` -> `(ndof,npoin)` transpose. Measured
   upside is only 1 to 2% in the realistic regime and it is an invasive,
   correctness-risky change. If the scatter ever becomes the clear
   remaining bottleneck after (1), it can be revisited, and even then the
   scatter-only 1.13x is the ceiling.
3. The RCM ordering is already capturing the available locality, so node
   renumbering is not a promising direction.

## Reproducing

The microbenchmark (`gs_bench.f90`) reads the real `ibool_sem2d.dat` and
times both layouts in both element orders. The `v` gather measurement uses
the `OPT_FINT_PROFILE` build flag with and without the gather guarded out.
