# An Heuristic Solver for the Mondrian Art Problem

### Introduction

The Mondrian Art Problem asks for a partition of the $n \times n$ square into multiple rectangles with integer side lengths, such that no two rectangles are congruent to each other (so if a $2 \times 3$ rectangle is already used, one can't place another $2 \times 3$ or a $3 \times 2$ rectangle). The challenge is to find for a given $n$ the smallest possible defect $M(n)$, which is defined as the difference between the area of the largest rectangle and smallest rectangle. You can try the problem yourself on my [interactive Mondrian Art editor](https://phoenixsmaug.github.io/mondrian-web/).

The optimal values are known up to $n = 66$, and are listed on OEIS as sequence [A276523](https://oeis.org/A276523). In 2016, Ed Pegg Jr published heuristic solutions up to $n = 900$ [on Wolfram Demonstrations](https://demonstrations.wolfram.com/MondrianArtProblem/) (technically his solution for $n = 196$ is invalid, since it contains the rectangle $26 \times 98$ twice). As far as I know, the source code used to find the solution is not public.

In this repository we present a Julia solver, which can achieve significantly lower defect solutions across the board:

![Record Comparison](img/record-comparison.png)

It also allows easily to extend heuristic solutions to $n = 3000$ and beyond:

![Range Extension](img/range-extension.png)

As an example, here is heuristic solution it finds for $n = 2000$:

![n = 2000 Example](img/mondrian_2000.png)

### How it works

Based on some public notes ([1](https://www.bradyharanblog.com/blog/mondrian-art-puzzle), [2](https://demonstrations.wolfram.com/MrsPerkinssQuilts/)), Ed Pegg Jr seems to have started with a explicit list of diverse rectangle partitions. In `src/generators.jl` we use the optimal solutions for small $n$ from OEIS (`data/import/oeis_solutions.txt`), Ed Pegg Jr's own solutions (`data/import/ed_pegg_jr.json`), random slicing trees based on [Guillotine cuts](https://en.wikipedia.org/wiki/Guillotine_partition) and electrical networks, whose suprising connection to packing problems is explained in [this YouTube video](https://youtu.be/0fH80JF2mDM). Note that we only save the topology, so which rectangles are adjacent to which other, and not the specific side lengths.

We can then apply a given rectangle partition to the $[0, 1] \times [0, 1]$ square and try to find a solution where all rectangles have the same area, but their side lengths are allowed to be real numbers. As explained in [this paper](https://arxiv.org/pdf/2007.09643), this transformation is just a big system of equations, which can easily be solved numerically with the Newton–Raphson method (see `src/geometry.jl`).

For a given $n$, we can then simply try to apply all of the cached $[0, 1] \times [0, 1]$ partitions by scaling them up and then rounding their internal wall coordinates to integers. The search space given by rounding up or down is then explored locally in `src/csp.jl`.

### Prime rectangulations (OEIS A340984)

Slicing trees, quilt substitutions and (to a lesser extent) electrical duals all tend to produce partitions that *factor*: some sub-collection of the rectangles has a rectangular union, so the partition is really a smaller partition with a block substituted into it. A partition where that never happens — no sub-collection of $2 \ldots k-1$ rectangles tiles a rectangle, so in particular no two rectangles share a complete face — is called **prime**, and those are counted by [A340984](https://oeis.org/A340984): $1, 1, 0, 0, 1, 0, 2, 6, 29, 119, 600, 3050, 16095, 87549$ cross-free classes for $k = 1 \ldots 14$. Prime partitions are exactly the ones the other generators cannot reach by construction, which makes them the natural missing seed family, and they matter for the Mondrian problem specifically: a shared face means two rectangles whose side lengths are locked together, which is precisely the rigidity a low-defect solution wants to avoid.

`data/import/a340984/` holds the complete catalogue of every class up to $k = 14$, taken from the [OEIS_A340984](https://github.com/PhoenixSmaug/OEIS_A340984) representative files. `src/io.jl` ingests them (`parse_a340984`), anisotropically normalizing each integer drawing onto $[0, 1]^2$; `src/prime.jl` can test any floorplan for primality (`is_prime_rectangulation`) and grow new ones.

Above $k = 14$ there is no catalogue, so `src/prime.jl` generates them. Splitting one rectangle of a prime partition in two always creates a block, so the split is followed by a short random walk of T-flips, re-equalizing after each move, until the block is gone; roughly 6% of attempts land on a prime partition with one more rectangle. Feeding each generated level back in as the parent pool for the next climbs the ladder as far as you like:

```bash
julia -t auto mondrian.jl gen-prime --kmax 20 --per-k 25000
```

On 12 threads that takes about 9 minutes and produces **253,446 prime classes** — the 107,452 complete ones through $k = 14$ plus 145,994 generated ones for $k = 15 \ldots 20$, de-duplicated up to the dihedral group of the square:

| $k$ | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|
| classes | 20,994 | 25,000 | 25,000 | 25,000 | 25,000 | 25,000 |
| rate | 335/s | 337/s | 309/s | 260/s | 273/s | 218/s |

The generated levels are a diverse sample, not an exhaustive enumeration — $a(15)$ alone is around half a million. The result is saved to `data/import/prime_generated.cache` and picked up automatically by `init-cache`, which uses the prime partitions in three places: as seeds for their own T-flip/wall-slide expansion, as substitution blocks for the quilt engine (the pinwheel and its relatives can *only* enter a larger partition by substitution), and as outer root templates for the hierarchical search.

Adding `--jsonl <dir>` also writes the generated partitions out in the catalogue's own schema, one `generated_nK.jsonl` per piece count, snapped onto the integer grid of their wall coordinates (`data/import/prime_generated_jsonl/`). Those files can be drawn with the catalogue's `render.py` exactly like the original representatives.

The crucial new idea over Ed's previous method is to use a two-level hierarchy of rectangles. So if you have a look at the solution for $n = 2000$ above, you'll notice that it is actually a partition of just $13$ big rectangles, each of which is split into $13$ smaller rectangles. As implemented in `src/solver.jl`, this two step optimization makes a search even with a lot of rectangles feasible. So while Ed's solution never used more than 27 rectangles, larger $n$ seems to benefit immensely from this increased number of rectangles.


### How to use the solver

Before running the solver for the first time, we need to built the cache of suitable rectangle partitions. With default options it takes around 45s of runtime and results in 115MB of data.
```bash
julia mondrian.jl init-cache
```

That uses the complete A340984 catalogue shipped in `data/import/a340984/`. To also seed it with prime partitions above $k = 14$, run the generator first — `init-cache` then finds `data/import/prime_generated.cache` on its own:
```bash
julia -t auto mondrian.jl gen-prime
```

Useful `init-cache` switches: `--prime-per-k <NUM>` caps how many prime classes are taken per piece count (`0` for all of them), `--expand-prime <NUM>` sets the expansion target per piece count, and `--no-prime` reproduces the old cache without them.

Afterwards, finding a heuristic solution for a given $n$, for example $213$ is just:
```bash
julia mondrian.jl 213
```

You can also supply `--deep` to extend the search caps and find better solutions:
```bash
julia mondrian.jl 213 --deep
```

Note that the specific defect of the solution found can vary by a few tiles, depending on the exact caps and cache initialization seeds used. The solver collects the best solution for each $n$ in `solutions.json`.

Using the short independent `verify.jl`, we can validate that the found solutions are valid solutions:
```
julia verify.jl solutions.json
```

Finally, you can use `plot.jl` to produce images of the solutions; like the $n = 2000$ solution seen above. Keep in mind that it needs further dependencies (`julia -e 'using Pkg; Pkg.add(["Luxor", "PicoSAT"])`).

### AI Declaration

This README is human-written, and all main ideas for the solver were proposed by humans. GPT 5.6 Sol helped optimizing the concrete implementation, especially the local CSP search and the electrical network generator. Additionally, Gemini 3.8 Flash was used to refactor and clean up code.
