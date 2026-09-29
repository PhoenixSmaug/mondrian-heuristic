/* enum_prime.c — exhaustive enumeration of cross-free prime rectangulations
 *
 * A prime rectangulation (OEIS A340984) is one in which no sub-collection of
 * 2..k-1 rectangles has a rectangular union.  Every cross-free class admits a
 * "tight" grid drawing in which each of its k-1 maximal wall segments sits on
 * its own grid line, so W + H = k + 1; enumerating those grids is complete.
 * Classes are identified by the catalogue's D4 wall-adjacency certificate,
 * because one class can have several distinct minimal drawings.
 *
 * Build: gcc -O3 -march=native -fopenmp -o enum_prime enum_prime.c
 * Usage: enum_prime <k> [threads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXK   24
#define MAXDIM (MAXK + 2)
#define MAXCELL (MAXDIM * MAXDIM)

typedef struct { int x0, y0, x1, y1; } Rect;

/* ---------- 128-bit hash set of certificates ---------------------------- */

typedef struct { uint64_t hi, lo; } Key;

typedef struct {
    Key     *slot;
    uint64_t mask;
    uint64_t count;
} Set;

static void set_init(Set *s, uint64_t cap_log2) {
    uint64_t n = 1ULL << cap_log2;
    s->slot = calloc(n, sizeof(Key));
    if (!s->slot) { fprintf(stderr, "out of memory\n"); exit(1); }
    s->mask = n - 1;
    s->count = 0;
}

static void set_grow(Set *s) {
    uint64_t oldn = s->mask + 1;
    Key *old = s->slot;
    uint64_t n = oldn * 2;
    s->slot = calloc(n, sizeof(Key));
    if (!s->slot) { fprintf(stderr, "out of memory on grow\n"); exit(1); }
    s->mask = n - 1;
    for (uint64_t i = 0; i < oldn; i++) {
        if (old[i].hi | old[i].lo) {
            uint64_t p = old[i].hi & s->mask;
            while (s->slot[p].hi | s->slot[p].lo) p = (p + 1) & s->mask;
            s->slot[p] = old[i];
        }
    }
    free(old);
}

/* returns 1 if newly inserted */
static int set_add(Set *s, Key k) {
    if (s->count * 10 >= (s->mask + 1) * 7) set_grow(s);
    uint64_t p = k.hi & s->mask;
    for (;;) {
        Key c = s->slot[p];
        if (!(c.hi | c.lo)) { s->slot[p] = k; s->count++; return 1; }
        if (c.hi == k.hi && c.lo == k.lo) return 0;
        p = (p + 1) & s->mask;
    }
}

/* ---------- certificate ------------------------------------------------- */

/* FNV-1a over a byte buffer, two independent seeds -> 128-bit key */
static Key hash128(const unsigned char *b, int n) {
    uint64_t h1 = 1469598103934665603ULL, h2 = 1099511628211ULL;
    for (int i = 0; i < n; i++) {
        h1 = (h1 ^ b[i]) * 1099511628211ULL;
        h2 = (h2 ^ b[i]) * 1469598103934665603ULL;
    }
    Key k = { h1 | 1, h2 };
    return k;
}

typedef struct {
    int side[4][MAXK];   /* right, bottom, left, top */
    int n[4];
} Nb;

/* BFS certificate of one oriented drawing, serialized into buf */
static int certificate(const Rect *r, int k, unsigned char *buf) {
    static _Thread_local Nb nb[MAXK];
    for (int i = 0; i < k; i++) {
        for (int s = 0; s < 4; s++) nb[i].n[s] = 0;
        /* right: s.x0 == r.x1, y-overlap; ordered by y0 */
        for (int j = 0; j < k; j++) {
            if (r[j].x0 == r[i].x1 && r[j].y0 < r[i].y1 && r[j].y1 > r[i].y0)
                nb[i].side[0][nb[i].n[0]++] = j;
        }
        for (int j = 0; j < k; j++) {
            if (r[j].y0 == r[i].y1 && r[j].x0 < r[i].x1 && r[j].x1 > r[i].x0)
                nb[i].side[1][nb[i].n[1]++] = j;
        }
        for (int j = 0; j < k; j++) {
            if (r[j].x1 == r[i].x0 && r[j].y0 < r[i].y1 && r[j].y1 > r[i].y0)
                nb[i].side[2][nb[i].n[2]++] = j;
        }
        for (int j = 0; j < k; j++) {
            if (r[j].y1 == r[i].y0 && r[j].x0 < r[i].x1 && r[j].x1 > r[i].x0)
                nb[i].side[3][nb[i].n[3]++] = j;
        }
        /* neighbours are collected in raster order of the rectangle list, which
         * is already sorted by (x0,y0); sort each side by its along-edge key */
        for (int s = 0; s < 4; s++) {
            int *a = nb[i].side[s], m = nb[i].n[s];
            for (int p = 1; p < m; p++) {
                int v = a[p], q = p - 1;
                int kv = (s == 0 || s == 2) ? r[v].y0 : r[v].x0;
                while (q >= 0 && ((s == 0 || s == 2) ? r[a[q]].y0 : r[a[q]].x0) > kv) {
                    a[q + 1] = a[q]; q--;
                }
                a[q + 1] = v;
            }
        }
    }

    int id[MAXK], order[MAXK];
    for (int i = 0; i < k; i++) id[i] = -1;
    int tl = -1;
    for (int i = 0; i < k; i++) if (r[i].x0 == 0 && r[i].y0 == 0) { tl = i; break; }
    if (tl < 0) return -1;
    id[tl] = 0; order[0] = tl;
    int cnt = 1;
    for (int qi = 0; qi < cnt; qi++) {
        int cur = order[qi];
        for (int s = 0; s < 4; s++)
            for (int p = 0; p < nb[cur].n[s]; p++) {
                int j = nb[cur].side[s][p];
                if (id[j] < 0) { id[j] = cnt; order[cnt++] = j; }
            }
    }
    if (cnt != k) return -1;

    int len = 0;
    for (int o = 0; o < k; o++) {
        int i = order[o];
        for (int s = 0; s < 4; s++) {
            buf[len++] = (unsigned char)nb[i].n[s];
            for (int p = 0; p < nb[i].n[s]; p++)
                buf[len++] = (unsigned char)id[nb[i].side[s][p]];
        }
    }
    return len;
}

/* the eight dihedral images of a drawing */
static void transform(const Rect *in, int k, int W, int H, int t, Rect *out) {
    for (int i = 0; i < k; i++) {
        int x0 = in[i].x0, y0 = in[i].y0, x1 = in[i].x1, y1 = in[i].y1;
        if (t & 4) { int nx0 = W - x1, nx1 = W - x0; x0 = nx0; x1 = nx1; }
        switch (t & 3) {
            case 0: out[i] = (Rect){x0, y0, x1, y1}; break;
            case 1: out[i] = (Rect){H - y1, x0, H - y0, x1}; break;
            case 2: out[i] = (Rect){W - x1, H - y1, W - x0, H - y0}; break;
            default: out[i] = (Rect){y0, W - x1, y1, W - x0}; break;
        }
    }
}

static int rect_cmp(const void *a, const void *b) {
    const Rect *p = a, *q = b;
    if (p->x0 != q->x0) return p->x0 - q->x0;
    if (p->y0 != q->y0) return p->y0 - q->y0;
    if (p->x1 != q->x1) return p->x1 - q->x1;
    return p->y1 - q->y1;
}

static Key canonical_key(const Rect *r, int k, int W, int H) {
    unsigned char best[8 * MAXK * 6], buf[8 * MAXK * 6];
    int bestlen = -1;
    Rect tmp[MAXK];
    for (int t = 0; t < 8; t++) {
        transform(r, k, W, H, t, tmp);
        qsort(tmp, k, sizeof(Rect), rect_cmp);
        int len = certificate(tmp, k, buf);
        if (len < 0) continue;
        if (bestlen < 0 || len < bestlen ||
            (len == bestlen && memcmp(buf, best, len) < 0)) {
            memcpy(best, buf, len); bestlen = len;
        }
    }
    if (bestlen < 0) { Key z = {0, 0}; return z; }
    return hash128(best, bestlen);
}

/* ---------- primality --------------------------------------------------- */

static int has_block(const Rect *r, int k, int W, int H) {
    /* A block is 2..k-1 pieces, so none exists below k = 3. Without this
     * guard the two halves of the k = 2 split share a full face and are
     * reported as a block, giving a(2) = 0 instead of 1. */
    if (k < 3) return 0;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < i; j++) {
            const Rect *a = &r[i], *b = &r[j];
            if (a->x0 == b->x0 && a->x1 == b->x1 && (a->y1 == b->y0 || b->y1 == a->y0)) return 1;
            if (a->y0 == b->y0 && a->y1 == b->y1 && (a->x1 == b->x0 || b->x1 == a->x0)) return 1;
        }
    /* exhaustive: every axis-aligned box that is exactly filled by 2..k-1 pieces */
    int xs[MAXDIM], ys[MAXDIM], nx = 0, ny = 0;
    char sx[MAXDIM] = {0}, sy[MAXDIM] = {0};
    for (int i = 0; i < k; i++) { sx[r[i].x0] = sx[r[i].x1] = 1; sy[r[i].y0] = sy[r[i].y1] = 1; }
    for (int i = 0; i <= W; i++) if (sx[i]) xs[nx++] = i;
    for (int i = 0; i <= H; i++) if (sy[i]) ys[ny++] = i;

    for (int a = 0; a < nx; a++)
        for (int b = a + 1; b < nx; b++)
            for (int c = 0; c < ny; c++)
                for (int d = c + 1; d < ny; d++) {
                    int x0 = xs[a], x1 = xs[b], y0 = ys[c], y1 = ys[d];
                    if (x0 == 0 && y0 == 0 && x1 == W && y1 == H) continue;
                    int area = 0, cnt = 0;
                    for (int i = 0; i < k; i++)
                        if (r[i].x0 >= x0 && r[i].x1 <= x1 && r[i].y0 >= y0 && r[i].y1 <= y1) {
                            cnt++; area += (r[i].x1 - r[i].x0) * (r[i].y1 - r[i].y0);
                        }
                    if (cnt >= 2 && cnt < k && area == (x1 - x0) * (y1 - y0)) return 1;
                }
    return 0;
}

/* ---------- raster fill ------------------------------------------------- */

typedef struct {
    int k, W, H;
    Rect r[MAXK];
    char cover[MAXCELL];
    Set *set;
    uint64_t leaves;
} Ctx;

static int shares_face(const Rect *a, const Rect *b) {
    if (a->x0 == b->x0 && a->x1 == b->x1 && (a->y1 == b->y0 || b->y1 == a->y0)) return 1;
    if (a->y0 == b->y0 && a->y1 == b->y1 && (a->x1 == b->x0 || b->x1 == a->x0)) return 1;
    return 0;
}

static void rec(Ctx *c, int placed, int filled) {
    int W = c->W, H = c->H, k = c->k;
    if (filled == W * H) {
        if (placed != k) return;
        /* minimality: every internal grid line must carry a rectangle edge */
        char sx[MAXDIM] = {0}, sy[MAXDIM] = {0};
        for (int i = 0; i < k; i++) { sx[c->r[i].x0] = sx[c->r[i].x1] = 1; sy[c->r[i].y0] = sy[c->r[i].y1] = 1; }
        for (int i = 1; i < W; i++) if (!sx[i]) return;
        for (int i = 1; i < H; i++) if (!sy[i]) return;
        if (has_block(c->r, k, W, H)) return;
        c->leaves++;
        Key key = canonical_key(c->r, k, W, H);
        if (key.hi | key.lo) set_add(c->set, key);
        return;
    }
    if (placed >= k) return;
    if (W * H - filled < k - placed) return;

    int fx = -1, fy = -1;
    for (int y = 0; y < H && fx < 0; y++)
        for (int x = 0; x < W; x++)
            if (!c->cover[y * W + x]) { fx = x; fy = y; break; }

    int maxw = 0;
    while (fx + maxw < W && !c->cover[fy * W + fx + maxw]) maxw++;

    for (int w = 1; w <= maxw; w++) {
        for (int h = 1; fy + h <= H; h++) {
            int ok = 1;
            for (int dx = 0; dx < w; dx++)
                if (c->cover[(fy + h - 1) * W + fx + dx]) { ok = 0; break; }
            if (!ok) break;

            Rect nr = { fx, fy, fx + w, fy + h };
            int bad = 0;
            if (k > 2)
                for (int i = 0; i < placed; i++)
                    if (shares_face(&nr, &c->r[i])) { bad = 1; break; }
            if (bad) continue;

            for (int dy = 0; dy < h; dy++)
                for (int dx = 0; dx < w; dx++) c->cover[(fy + dy) * W + fx + dx] = 1;
            c->r[placed] = nr;
            rec(c, placed + 1, filled + w * h);
            for (int dy = 0; dy < h; dy++)
                for (int dx = 0; dx < w; dx++) c->cover[(fy + dy) * W + fx + dx] = 0;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <k> [threads]\n", argv[0]); return 1; }
    int k = atoi(argv[1]);
    if (k < 1 || k > MAXK) { fprintf(stderr, "k out of range\n"); return 1; }
#ifdef _OPENMP
    if (argc > 2) omp_set_num_threads(atoi(argv[2]));
#endif
    double t0 = (double)clock() / CLOCKS_PER_SEC;
    struct timespec ts0, ts1;
    timespec_get(&ts0, TIME_UTC);

    /* Tasks are (tight grid, first rectangle at the top-left cell). One grid
     * dominates the work, so splitting only by grid leaves 11 threads idle. */
    int task[MAXDIM * MAXDIM * MAXDIM][4], nt = 0;
    for (int W = 1; W <= k; W++) {
        int H = k + 1 - W;
        if (H < W || H < 1) continue;
        if (W * H < k) continue;
        for (int w = 1; w <= W; w++)
            for (int h = 1; h <= H; h++) {
                if (W * H - w * h < k - 1) continue;
                task[nt][0] = W; task[nt][1] = H; task[nt][2] = w; task[nt][3] = h; nt++;
            }
    }

    Set global; set_init(&global, 20);
    uint64_t leaves = 0;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int p = 0; p < nt; p++) {
        Ctx c;
        memset(&c, 0, sizeof c);
        c.k = k; c.W = task[p][0]; c.H = task[p][1];
        int w = task[p][2], h = task[p][3];
        Set local; set_init(&local, 14);
        c.set = &local;
        c.r[0] = (Rect){0, 0, w, h};
        for (int dy = 0; dy < h; dy++)
            for (int dx = 0; dx < w; dx++) c.cover[dy * c.W + dx] = 1;
        rec(&c, 1, w * h);
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            for (uint64_t i = 0; i <= local.mask; i++)
                if (local.slot[i].hi | local.slot[i].lo) set_add(&global, local.slot[i]);
            leaves += c.leaves;
        }
        free(local.slot);
    }

    timespec_get(&ts1, TIME_UTC);
    double wall = (ts1.tv_sec - ts0.tv_sec) + 1e-9 * (ts1.tv_nsec - ts0.tv_nsec);
    double cpu = (double)clock() / CLOCKS_PER_SEC - t0;
    printf("k=%d  cross-free classes=%llu  leaves=%llu  wall=%.2fs cpu=%.2fs\n",
           k, (unsigned long long)global.count, (unsigned long long)leaves, wall, cpu);
    return 0;
}
