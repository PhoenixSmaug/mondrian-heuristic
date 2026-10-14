# An Heuristic Solver for the Mondrian Art Problem

### Introduction

The Mondrian Art Problem asks us to cut an $n \times n$ square into rectangles with integer side lengths so that no two rectangles are congruent: once a $2 \times 3$ rectangle is used, neither another $2 \times 3$ nor a $3 \times 2$ may appear. The *defect* of such a partition is the area of its largest rectangle minus the area of its smallest, and the goal is the smallest defect possible for a given $n$, written $M(n)$. You can try the problem yourself on the [interactive Mondrian Art editor](https://phoenixsmaug.github.io/mondrian-web/).

Optimal values are known only up to $n = 66$, and are listed on OEIS as [A276523](https://oeis.org/A276523). For larger $n$ we can only find *heuristic solutions*: valid partitions found by search, with no proof that a smaller defect is impossible. A heuristic solution therefore gives an upper bound on $M(n)$, not its exact value. In 2016 Ed Pegg Jr published heuristic solutions up to $n = 900$ [on Wolfram Demonstrations](https://demonstrations.wolfram.com/MondrianArtProblem/) (his solution for $n = 196$ is in fact invalid, since it uses the rectangle $26 \times 98$ twice); as far as we know, the code behind them is not public.

This repository presents a Julia solver that finds lower defects than Ed's across the board and extends heuristic solutions to $n = 3000$ and beyond:

![Ed Pegg Jr and this solver](img/census.png)

The values in this graph and their partitions, each re-verified by an independent checker, are in `records/`: `values.txt` lists the best defect for every $n$ up to 3000, `solutions.json` holds the partitions found by the current version, and `census.jsonl` gives, for every $n$ from 67 to 3000, the defect found by each of its two schedules. Those runs used 240 s per $n$, at which the timeout cuts the search short for $n$ above roughly 1500, so the large-$n$ values understate what the solver finds with more time.

As an example, here is the solution it finds for $n = 2000$:

![n = 2000 Example](img/mondrian_2000.png)

### How it works

Judging by public notes ([1](https://www.bradyharanblog.com/blog/mondrian-art-puzzle), [2](https://demonstrations.wolfram.com/MrsPerkinssQuilts/)), Ed Pegg Jr started from an explicit list of diverse rectangle partitions, and this solver does the same. The list is built in `src/generators.jl` from the optimal solutions for small $n$ (`data/import/oeis_solutions.txt`), Ed Pegg Jr's own solutions (`data/import/ed_pegg_jr.json`), random slicing trees based on [guillotine cuts](https://en.wikipedia.org/wiki/Guillotine_partition), and electrical networks, whose surprising connection to packing problems is explained in [this video](https://youtu.be/0fH80JF2mDM). Only the topology of each partition is stored, meaning which rectangle is adjacent to which, not its side lengths.

A topology is placed on the unit square, and the solver looks for real side lengths that give every rectangle the same area. As explained in [this paper](https://arxiv.org/pdf/2007.09643), that is a system of equations, solved numerically with the Newton–Raphson method in `src/geometry.jl`. For a given $n$ the solver scales these partitions up and rounds their interior walls to integers; the choices of rounding each wall up or down are then explored by a local constraint search, the CSP in `src/csp.jl`.

The crucial idea, beyond Ed's method, is a two-level hierarchy. The $n = 2000$ solution above is a partition into just $13$ big rectangles, each of which is split into $13$ smaller ones. Optimizing the two levels separately makes a search with many rectangles feasible. Ed's solutions never used more than 27 rectangles, and larger $n$ benefit greatly from using more.

### Prime rectangulations (OEIS A340984)

Slicing trees, quilt substitutions and, to a lesser extent, electrical networks tend to produce partitions that *factor*: some sub-collection of the rectangles has a rectangular union, so the partition is really a smaller partition with a block substituted into it. A partition of $k$ rectangles in which no $2$ to $k-1$ of them tile a rectangle is called *prime*. The smallest example is the pinwheel, with five rectangles around a central one. Prime partitions are counted by [A340984](https://oeis.org/A340984): $1, 1, 0, 0, 1, 0, 2, 6, 29, 119, 600, 3054, 16143, 88112$ for $k = 1$ to $14$. From $k = 12$ on, a few of these are *crossed*, with four rectangles meeting at a point; the solver uses only the cross-free ones, of which there are $3050$, $16095$ and $87549$ for $k = 12, 13, 14$.

Prime partitions are exactly the ones the other generators cannot produce by construction, which makes them the natural missing family of seeds. They also matter for this problem in particular: two rectangles that share a complete side have their side lengths locked together, and that is precisely the rigidity a low-defect solution needs to avoid.

The complete catalogue up to $k = 14$, the files `representatives_n1.jsonl` to `representatives_n14.jsonl` from [meghanto/OEIS_A340984](https://github.com/meghanto/OEIS_A340984), is included compressed as `data/import/a340984.tar.xz` (3 MB instead of 81 MB). The function `parse_a340984` in `src/io.jl` reads them and normalizes each drawing onto the unit square, and `src/prime.jl` can test any partition for primality (`is_prime_rectangulation`) and grow new ones.

Above $k = 14$ there is no catalogue, so the solver generates prime partitions itself. Splitting one rectangle of a prime partition in two always creates a block, so the split is followed by a short random walk of T-flips, re-equalizing after each move, until the block is gone; about 6% of attempts give a prime partition with one more rectangle. Each generated level becomes the parent pool for the next:

```bash
julia -t auto mondrian.jl gen-prime --kmax 20 --per-k 25000
```

On 12 threads this takes about 9 minutes and gives 253,446 prime classes, counted up to rotation and reflection: the 107,452 cross-free classes through $k = 14$, and 145,994 generated ones for $k = 15$ to $20$.

| $k$ | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|
| classes | 20,994 | 25,000 | 25,000 | 25,000 | 25,000 | 25,000 |

The generated levels are a diverse sample, not a complete enumeration; the count for $k = 15$ alone is around half a million. They are saved to `data/import/prime_generated.cache`, which the cache builder picks up automatically and uses in three ways: as seeds for its own expansion, as blocks for quilt substitution (the pinwheel and its relatives can enter a larger partition only by substitution), and as outer templates for the two-level search. Adding `--jsonl <dir>` also writes them in the catalogue's own format, so the catalogue's `render.py` can draw them.

### An evolved decision layer

The machinery above produces candidate partitions; the decision layer in `src/solver.jl` decides how to spend effort on them. It chooses which hierarchy shapes to try, where a shape $(R, M)$ means $R$ outer rectangles each split into $M$ inner ones; which outer roots and rounding masks get CSP calls; and how a fixed budget of CSP calls is shared out. We evolved this layer with [ShinkaEvolve](https://github.com/SakanaAI/ShinkaEvolve), an LLM-driven program-evolution loop, and left everything below it untouched. Each candidate program was scored on 30 values of $n$ between 300 and 3000, at 900 s per $n$, by the quantity
$$50 \times \text{mean}\left(\frac{d_{\text{best known}}}{d_{\text{found}}}\right),$$
so that a score of 50 means matching the best-known values everywhere. The winners were then scored again on 50 other values of $n$ they had never been tuned on.

**Successive halving.** The budget for each $n$ is 100 CSP calls, and there are about 16 candidate shapes, called *arms*. Successive halving gives every arm a small equal share of calls, ranks the arms by the best defect each has found, keeps the better half and repeats with a larger share for each survivor. Arms whose share is not enough for a valid partition yet are ranked by an estimated defect that costs no calls. The arms are *resumable*: an arm that survives a round continues from where it stopped, so nothing it explored is repeated. In the default schedule, for example, 16 arms get 2 calls each in the first round, the 8 survivors 4 each, and the 4 after that 8 each, about 32 calls per round, so most of the budget ends up on the shapes that are actually producing low defects.

Two evolved programs are kept because they win on different $n$. They share all their machinery and differ only in when the scheduler adds new arms, namely the neighbouring shapes $(R \pm 1, M)$ and $(R, M \pm 1)$ of the best arm:

| `--schedule` | when neighbouring arms are added | evolution set | held-out 50 |
|---|---|---|---|
| `bandit` (default) | as soon as two or fewer arms survive; they join the next round of halving | 55.95 | 55.23, best known beaten on 47 of 50 |
| `race` | only after the pool has halved to one survivor and every arm has run out of roots | 55.50 | not run |
| `best` | runs both schedules and keeps the lower defect, taking twice the time | | |

Evolution started from our earlier rewrite of `solver.jl`, which scored 52.42 on the same held-out values (39 of 50 beaten). That rewrite already chose shapes near a fitted piece-count law of about $2.895\,n^{0.472}$ rectangles, ranked outer roots by how evenly they split the square before spending any call on them, and counted CSP calls as its budget. Evolution added the resumable arms, successive halving, the widening to neighbouring shapes, and an ordering of shapes by their measured quality in an earlier survey. The evolved code also contains a non-uniform arm that never runs, because it calls `equalize_weighted` with floating-point weights, for which there is no method; we leave it in so that the file reproduces the evaluated programs exactly.

Because the budget is a count of calls, a result is deterministic as long as `--timeout` does not cut the search short. At 240 s per $n$ the timeout does cut it short for $n$ above roughly 1500, where results then depend on the speed of the machine. Even at 900 s it can at the largest $n$: with `race`, the value $n = 2713$ gave 265 in one run and 280 in another.

### How to use the solver

Before the first run, unpack the A340984 catalogue and build the cache of rectangle partitions:
```bash
tar -xJf data/import/a340984.tar.xz -C data/import
```

Building the cache with the default options takes about 45 s and produces 115 MB of data. To include the generated prime partitions above $k = 14$, run the generator first; the cache builder then finds its output on its own:
```bash
julia -t auto mondrian.jl gen-prime
julia mondrian.jl init-cache
```

Useful options for the cache builder are `--prime-per-k <NUM>`, which caps the prime classes taken per piece count (`0` for all of them), `--expand-prime <NUM>`, which sets the expansion target per piece count, and `--no-prime`, which leaves the primes out. The generator is multi-threaded and time-limited, so its output, and with it the cache, varies from machine to machine. To use exactly the cache behind `records/`, download `floorplans.cache` and `outer_roots.cache` from this repository's release assets into `data/cache/` instead of building one.

Finding a heuristic solution for a given $n$, for example $213$, is then just:
```bash
julia mondrian.jl 213
```

To run both evolved schedules and keep the better result:
```bash
julia mondrian.jl 213 --schedule best
```

The flag `--deep` raises the budget from 100 to 3600 CSP calls; we have not yet evaluated it for the evolved layer:
```bash
julia mondrian.jl 213 --deep
```

The solver collects the best solution found for each $n$ in `solutions.json`. The short, independent `verify.jl` checks that they are valid:
```bash
julia verify.jl solutions.json
```

Finally, `plot.jl` draws solutions such as the $n = 2000$ example above. It needs two further packages, installed with `julia -e 'using Pkg; Pkg.add(["Luxor", "PicoSAT"])'`.

### AI Declaration

The original solver and the original version of this README were written by Mia Muessig, and all main ideas for that solver were proposed by humans. GPT 5.6 Sol helped optimize the implementation, especially the local CSP search and the electrical network generator, and Gemini 3.8 Flash was used to refactor and clean up code.

The later additions were directed by Thomas Garrison, who proposed many of their ideas, and were written with Claude, partly inside the ShinkaEvolve loop; `enum/certificate.py` is copied from OEIS_A340984. This README, in its present form, was written by Claude.
