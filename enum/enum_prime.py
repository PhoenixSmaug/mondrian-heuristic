"""Exhaustive enumeration of prime (A340984) rectangulations with k rectangles.

Every class has a unique minimal grid drawing: collapse the distinct wall
coordinates to 0..W and 0..H. Minimality forces every internal grid line to
carry at least one rectangle edge, which bounds W+H <= k+1. So enumerating all
W x H grids in that range, raster-filling them with exactly k rectangles, and
keeping the prime ones is complete.
"""
import sys
from functools import lru_cache

from certificate import canonical_d4_certificate

def canonical(rects, W, H):
    """Topology invariant, not a drawing: the same rectangulation can have
    several distinct minimal grid drawings, so dedup has to happen at the
    wall-adjacency level (the catalogue's own d4_certificate)."""
    return canonical_d4_certificate(list(rects), W, H)

def has_block(rects, k):
    xs = sorted({r[0] for r in rects} | {r[2] for r in rects})
    ys = sorted({r[1] for r in rects} | {r[3] for r in rects})
    for i1 in range(len(xs)):
        for j1 in range(len(ys)):
            for i2 in range(i1 + 1, len(xs)):
                for j2 in range(j1 + 1, len(ys)):
                    x0, y0, x1, y1 = xs[i1], ys[j1], xs[i2], ys[j2]
                    if (x1 - x0) * (y1 - y0) >= sum((r[2]-r[0])*(r[3]-r[1]) for r in rects):
                        continue
                    ids = [r for r in rects if x0 <= r[0] and r[2] <= x1 and y0 <= r[1] and r[3] <= y1]
                    if 2 <= len(ids) < k and sum((r[2]-r[0])*(r[3]-r[1]) for r in ids) == (x1-x0)*(y1-y0):
                        return True
    return False

def shares_full_face(a, b):
    if a[0] == b[0] and a[2] == b[2] and (a[3] == b[1] or b[3] == a[1]): return True
    if a[1] == b[1] and a[3] == b[3] and (a[2] == b[0] or b[2] == a[0]): return True
    return False

def fill(W, H, k, classes):
    cover = [[False]*W for _ in range(H)]
    rects = []
    total = W*H

    def first_free():
        for y in range(H):
            row = cover[y]
            for x in range(W):
                if not row[x]:
                    return x, y
        return None

    def used_lines():
        vx = {r[0] for r in rects} | {r[2] for r in rects}
        vy = {r[1] for r in rects} | {r[3] for r in rects}
        return all(i in vx for i in range(1, W)) and all(j in vy for j in range(1, H))

    def rec(placed, filled):
        if filled == total:
            if placed == k and used_lines() and not has_block(rects, k):
                classes.add(canonical(rects, W, H))
            return
        if placed >= k:
            return
        # every remaining piece covers >= 1 cell, and <= (k-placed) pieces remain
        if total - filled < k - placed:
            return
        x, y = first_free()
        maxw = 0
        while x + maxw < W and not cover[y][x + maxw]:
            maxw += 1
        for w in range(1, maxw + 1):
            h = 1
            while y + h <= H:
                if any(cover[y + h - 1][x + dx] for dx in range(w)):
                    break
                new = (x, y, x + w, y + h)
                if k > 2 and any(shares_full_face(new, r) for r in rects):
                    h += 1
                    continue
                for dy in range(h):
                    for dx in range(w):
                        cover[y + dy][x + dx] = True
                rects.append(new)
                rec(placed + 1, filled + w * h)
                rects.pop()
                for dy in range(h):
                    for dx in range(w):
                        cover[y + dy][x + dx] = False
                h += 1
    rec(0, 0)

def enumerate_k(k):
    classes = set()
    for W in range(1, k + 1):
        for H in range(1, k + 2 - W):
            if W * H < k:
                continue
            fill(W, H, k, classes)
    return classes

if __name__ == "__main__":
    import time
    for k in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        t0 = time.time()
        c = enumerate_k(k)
        print(f"k={k:3d}  classes={len(c):9d}   {time.time()-t0:8.2f}s", flush=True)
