#!/usr/bin/env python3
"""
A real, if small, 2-D surface-wave travel-time tomography, for resolution tests.

Why this exists
---------------
`synthetic.py` fakes a degraded inversion by mixing the true field with a blurred
copy, weighted by distance from the array centre.  That is fine for exercising
plumbing, but it is NOT a resolution test, and for judging a distance-weighted
merge it is circular: the imposed error structure has the same functional form as
the merge weights, so distance weighting is guaranteed to look good.

Here nothing about the degradation is imposed.  Synthetic travel times are
computed through a known model along real station-pair paths, noise is added, and
the model is recovered by damped, smoothed least squares on the same grid.
Amplitude loss, smearing along dominant path azimuths, edge artefacts and the
null space all emerge from the coverage and the regularisation, exactly as in a
published inversion (Aki & Lee, 1976; Barmin et al., 2001; Rawlinson & Spakman,
2016).  Whatever the merge then achieves, it achieves against errors it did not
help design.

numpy only.  The normal equations are formed directly, which also yields the
model resolution matrix needed for the standard diagnostics.
"""
import numpy as np

DEG_KM = 111.195


# ============================================================================
class Grid:
    """Regular lon/lat cell grid.  Cells are indexed [i_lon, j_lat] -> flat."""

    def __init__(self, lon0, lon1, lat0, lat1, d):
        self.d = float(d)
        self.lonv = np.round(np.arange(lon0, lon1 + 1e-9, d), 6)
        self.latv = np.round(np.arange(lat0, lat1 + 1e-9, d), 6)
        self.nx, self.ny = len(self.lonv), len(self.latv)
        self.n = self.nx*self.ny
        self.lat_mid = float(np.mean(self.latv))
        self.kx = np.cos(np.deg2rad(self.lat_mid))     # lon degrees -> equivalent

    def flat(self, i, j):
        return i*self.ny + j

    def index(self, lon, lat):
        """Flat cell index for arrays of coordinates; -1 outside the grid."""
        i = np.rint((np.asarray(lon) - self.lonv[0])/self.d).astype(int)
        j = np.rint((np.asarray(lat) - self.latv[0])/self.d).astype(int)
        bad = (i < 0) | (i >= self.nx) | (j < 0) | (j >= self.ny)
        out = i*self.ny + j
        out[bad] = -1
        return out

    def to_2d(self, m):
        return np.asarray(m, float).reshape(self.nx, self.ny)

    def xyz(self, m):
        X, Y = np.meshgrid(self.lonv, self.latv, indexing="ij")
        return X.ravel(), Y.ravel(), np.asarray(m, float).ravel()

    def laplacian_normal(self):
        """L^T L for a first-difference smoothness penalty between neighbours."""
        LtL = np.zeros((self.n, self.n))
        for i in range(self.nx):
            for j in range(self.ny):
                a = self.flat(i, j)
                for di, dj in ((1, 0), (0, 1)):
                    if i + di < self.nx and j + dj < self.ny:
                        b = self.flat(i + di, j + dj)
                        LtL[a, a] += 1.0; LtL[b, b] += 1.0
                        LtL[a, b] -= 1.0; LtL[b, a] -= 1.0
        return LtL


# ============================================================================
def stations(grid, n=170, seed=11, jitter=0.45):
    """A quasi-uniform station distribution: a jittered grid, like a real array."""
    rng = np.random.default_rng(seed)
    aspect = (grid.lonv[-1] - grid.lonv[0])*grid.kx / (grid.latv[-1] - grid.latv[0])
    ny = max(2, int(round(np.sqrt(n/aspect))))
    nx = max(2, int(round(n/ny)))
    lo = np.linspace(grid.lonv[0] + .6, grid.lonv[-1] - .6, nx)
    la = np.linspace(grid.latv[0] + .6, grid.latv[-1] - .6, ny)
    X, Y = np.meshgrid(lo, la, indexing="ij")
    st = np.c_[X.ravel(), Y.ravel()]
    st[:, 0] += rng.uniform(-jitter, jitter, len(st))/grid.kx
    st[:, 1] += rng.uniform(-jitter, jitter, len(st))
    return st


def station_pairs(st, sep_min=1.5, sep_max=10.0, kx=1.0):
    """Every station pair whose separation lies in the usable range."""
    out = []
    for a in range(len(st)):
        for b in range(a + 1, len(st)):
            dx = (st[b, 0] - st[a, 0])*kx
            dy = st[b, 1] - st[a, 1]
            s = np.hypot(dx, dy)
            if sep_min <= s <= sep_max:
                out.append((a, b, s))
    return out


def ray_rows(grid, st, pairs, step=0.08):
    """
    Sparse rows of the design matrix: for each path, the length IN KILOMETRES
    travelled inside each grid cell.  A straight path in the local scaled plane,
    sampled finely and binned - adequate at these scales and free of the corner
    cases of exact cell clipping.

    Kilometres matter: with slowness in s/km the data are travel times in
    seconds, so the noise level can be set to a real dispersion-measurement
    uncertainty rather than an arbitrary number.
    """
    rows, lengths, keys = [], [], []
    for a, b, s in pairs:
        nseg = max(2, int(np.ceil(s/step)))
        t = (np.arange(nseg) + 0.5)/nseg                 # segment midpoints
        lon = st[a, 0] + (st[b, 0] - st[a, 0])*t
        lat = st[a, 1] + (st[b, 1] - st[a, 1])*t
        idx = grid.index(lon, lat)
        idx = idx[idx >= 0]
        if len(idx) < 2:
            continue
        seg = s*DEG_KM/nseg              # kilometres
        u, c = np.unique(idx, return_counts=True)
        rows.append((u, c.astype(float)*seg))
        lengths.append(s)
        # identity of the measurement, not of the subarray that used it: two
        # overlapping subarrays share stations, so they share paths, and a shared
        # path is ONE measurement carrying ONE noise realisation
        keys.append(pair_key(st[a], st[b]))
    return rows, np.asarray(lengths), keys


def pair_key(pa, pb):
    """Order-independent identity of a station pair, to 1e-4 degrees."""
    a = (round(float(pa[0]), 4), round(float(pa[1]), 4))
    b = (round(float(pb[0]), 4), round(float(pb[1]), 4))
    return (a, b) if a <= b else (b, a)


class MeasurementNoise:
    """
    One noise realisation per station pair, shared between subarrays.

    Drawing noise per inversion would give the same physical measurement two
    different errors in the two subarrays that use it, which makes the pair of
    models look more independent than it is and flatters any averaging.
    """

    def __init__(self, sigma, seed=0):
        self.sigma = float(sigma)
        self.seed = int(seed)
        self._v = {}

    def vector(self, keys):
        out = np.empty(len(keys))
        for i, k in enumerate(keys):
            if k not in self._v:
                h = abs(hash((k, self.seed))) % (2**32)
                self._v[k] = float(np.random.default_rng(h).standard_normal())
            out[i] = self.sigma*self._v[k]
        return out


def forward(rows, slowness):
    """Travel times = integral of slowness along each path."""
    s = np.asarray(slowness, float)
    return np.array([float(v @ s[u]) for u, v in rows])


def normal_equations(rows, n):
    """G^T G and the operator needed to form G^T d, accumulated row by row."""
    GtG = np.zeros((n, n))
    for u, v in rows:
        GtG[np.ix_(u, u)] += np.outer(v, v)
    return GtG


def gtd(rows, d, n):
    out = np.zeros(n)
    for (u, v), di in zip(rows, d):
        out[u] += v*di
    return out


def ray_density(rows, n):
    """Total path length in each cell: the standard coverage diagnostic."""
    out = np.zeros(n)
    for u, v in rows:
        out[u] += v
    return out


# ============================================================================
class Inversion:
    """
    Damped, smoothed linear inversion on a fixed grid and path set.

    Minimises  ||G m - d||^2 + alpha^2 ||m - m0||^2 + beta^2 ||L m||^2,
    the standard regularised tomography formulation.  Factorising once lets the
    same operator invert many synthetic data sets cheaply, and gives the model
    resolution matrix R = (G^T G + alpha^2 I + beta^2 L^T L)^-1 G^T G.
    """

    def __init__(self, grid, rows, alpha=6.0, beta=40.0):
        self.grid, self.rows = grid, rows
        self.alpha, self.beta = float(alpha), float(beta)
        n = grid.n
        self.GtG = normal_equations(rows, n)
        A = self.GtG + (self.alpha**2)*np.eye(n) + (self.beta**2)*grid.laplacian_normal()
        self.Ainv = np.linalg.inv(A)
        self.density = ray_density(rows, n)

    def invert(self, d, m0):
        m0 = np.asarray(m0, float)
        rhs = gtd(self.rows, np.asarray(d, float) - forward(self.rows, m0), self.grid.n)
        return m0 + self.Ainv @ rhs

    def resolution(self):
        """Model resolution matrix.  Its diagonal is the usual resolution map."""
        return self.Ainv @ self.GtG

    def recover(self, true_slowness, m0, noise=0.0, seed=0, noise_vec=None,
                fwd_rows=None, fwd_slowness=None):
        """
        Forward-model, add noise, invert: one complete synthetic experiment.

        fwd_rows / fwd_slowness let the data be generated on a FINER grid than
        the one inverted for.  Computing the synthetic travel times through the
        same discretisation that is then solved for is an inverse crime: the true
        model is exactly representable in the model space and recovery comes out
        flattered.  Passing a finer forward grid removes that.

        noise_vec overrides the internal draw, so that a measurement shared by
        two subarrays carries the same error in both.
        """
        if fwd_rows is not None and fwd_slowness is not None:
            d = forward(fwd_rows, fwd_slowness)
        else:
            d = forward(self.rows, true_slowness)
        if noise_vec is not None:
            d = d + np.asarray(noise_vec, float)
        elif noise:
            d = d + noise*np.random.default_rng(seed).standard_normal(len(d))
        return self.invert(d, m0)

    def well_covered(self, frac=0.12):
        """
        Mask of cells a study would actually publish: those whose ray coverage
        exceeds a fraction of the median coverage.  Real papers mask their maps
        this way, and it makes the two subarray footprints ragged rather than
        rectangular, which is the realistic case for a merge.
        """
        med = np.median(self.density[self.density > 0]) if np.any(self.density > 0) else 0.0
        return self.density >= frac*med


# ============================================================================
def checkerboard(grid, cell=3.0, amp=0.12, base=3.85, pattern="blocks"):
    """A checkerboard velocity model, in km/s, on the grid."""
    X, Y = np.meshgrid(grid.lonv, grid.latv, indexing="ij")
    b = (np.sin(np.pi*(X - grid.lonv[0])/cell)
         * np.sin(np.pi*(Y - grid.latv[0])/cell))
    if pattern == "blocks":
        b = np.sign(b)
    return (base + amp*b).ravel()


def spike(grid, at, amp=0.12, base=3.85):
    """A single-cell anomaly, for a point-spread-function test."""
    m = np.full(grid.n, base)
    k = int(grid.index(np.array([at[0]]), np.array([at[1]]))[0])
    if k >= 0:
        m[k] += amp
    return m


def amplitude_recovery(true_v, rec_v, mask, base=None):
    """
    Least-squares slope of recovered anomaly on true anomaly, over `mask`.
    1.0 is full recovery; regularised inversions return less.
    """
    t = np.asarray(true_v, float)[mask]
    r = np.asarray(rec_v, float)[mask]
    b = t.mean() if base is None else base
    dt, dr = t - b, r - r.mean()
    return float((dt @ dr)/(dt @ dt)) if dt @ dt > 0 else np.nan
