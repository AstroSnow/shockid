# -*- coding: utf-8 -*-
"""
shockID: automated detection and classification of MHD shocks in 2D and 3D simulation data.

Pipeline (see shockid()):
  1. optionally smooth the data (smoothdata)
  2. candidate cells = converging flow, div(v) < -convl (divergence)
  3. thin candidates to one cell per front with non-maximum suppression of |grad rho|
     along the density-gradient normal (nonMaxSuppress)
  4. at each remaining cell sample the variables along the normal (getNormVals), optionally
     refine the normal with div(B)=0 (refineNormalDB), move to the shock frame
     (getShockFrame) and compare the normal flow speed with the slow / Alfven / fast speeds
     either side (classifyShock, getWaveSpeeds, getState)
  5. the pre->post state pair gives the shock type (TRANS)

Array conventions: 2D data are indexed [y, x], 3D data [z, y, x]. In the code col = y index,
row = x index and zrow = z index. Returned shock locations are [y, x] or [z, y, x].

Units: wave speeds use vA^2 = B^2/rho (Heaviside-Lorentz, mu0 = 1). For Gaussian data
(magnetic pressure B^2/8pi) pass B/sqrt(4*pi).
"""
import numpy as np
import time
from scipy.ndimage import gaussian_filter, map_coordinates

# (pre-shock state, post-shock state) -> shock type. States are defined in getState:
# 1 super-fast, 2 between Alfven and fast, 3 between slow and Alfven, 4 sub-slow.
TRANS = {(1,2):'fast', (3,4):'slow', (1,3):'int1',
         (1,4):'int2', (2,3):'int3', (2,4):'int4'}


###############################################################################
def shockid(gridx, gridy, gridz, rog, vxg, vyg, vzg, bxg, byg, bzg, prg, ndim=2, smthfac=0,
            nproc=1, convl=0.0001, avecyl=2, normsmooth=0.0, dbproject=False, periodic=True,
            shocktol=0.05, nbrsearch=1):
    """Find and classify shocks.

    gridx, gridy, gridz : 1D coordinate arrays (gridz ignored in 2D). Grids may be
                 non-uniform (stretched): all derivatives use the coordinates.
    rog, vxg, vyg, vzg, bxg, byg, bzg, prg : density, velocity, magnetic field and gas
                 pressure, indexed [y,x] (2D) or [z,y,x] (3D). In 2D, vzg and bzg may be the
                 scalar 0.0.
    ndim       : 2 or 3.
    smthfac    : Gaussian sigma (cells) applied to all variables before anything else. 0 = off.
    nproc      : number of processes for the classification loop.
    convl      : candidate cells need div(v) < -convl.
    avecyl     : half-length (cells) of the line sampled along the normal; the pre/post
                 states are taken at its ends.
    normsmooth : Gaussian sigma (cells) applied to density ONLY for computing the shock normal
                 (and the non-maximum suppression). Jump values still come from the
                 smthfac-smoothed data. 0 = off.
    dbproject  : refine each normal using div(B)=0: remove its component along
                 dB = B_post - B_pre and resample. Skipped where |dB| < 5% of |B|.
    periodic   : True  - sampling windows near the boundary wrap around (periodic domains).
                 False - discard candidates within max(avecyl,2) cells of any boundary,
                         so no window ever leaves the domain (non-periodic domains).
    shocktol   : any Mach number within shocktol of 1 makes that side 'marginal' (state 0)
                 and the candidate is dropped.
    nbrsearch  : cells along the normal to try when the gradient maximum itself cannot be
                 classified. The result is recorded at the maximum, so nothing is double
                 counted. 0 = off.

    Returns a dict with, for each type k in TRANS.values():
      k        : (N, ndim) int array of locations, [y,x] or [z,y,x]
      k+'mach' : (N, 6) slow, Alfven, fast Mach numbers, pre-shock then post-shock
      k+'qual' : (N, 4) [mach_margin, bn_jump, bt_rotation_deg, shift], see shockQuality
                 and shockClassLoop.
      k+'norm' : (N, ndim) unit shock normal used for the classification (after any
                 dbproject refinement), oriented from the pre- to the post-shock side.
      k+'perp' : (N, ndim) unit tangent to the front, see shockVectors. In 3D the second
                 tangent is np.cross(norm, perp).
                 Vector components are (x, y[, z]), the reverse of the [z,]y,x order of
                 the locations.
    """
    # Coordinates (not spacings), so derivatives are correct on stretched grids
    x = np.asarray(gridx)
    y = np.asarray(gridy)
    z = np.asarray(gridz) if ndim == 3 else 0.0

    ds = smoothdata(rog, vxg, vyg, vzg, bxg, byg, bzg, prg, smthfac)

    print('Calculating velocity divergence')
    divv = divergence(ds, ndim, x, y, z)

    # Density gradients define the shock normal
    print('Calculating density gradients')
    ronorm = gaussian_filter(ds['ro'], normsmooth) if normsmooth > 0 else ds['ro']
    gradroz = 0.0
    if ndim == 2:
        gradrox = np.gradient(ronorm, x, axis=1)
        gradroy = np.gradient(ronorm, y, axis=0)
    else:
        gradrox = np.gradient(ronorm, x, axis=2)
        gradroy = np.gradient(ronorm, y, axis=1)
        gradroz = np.gradient(ronorm, z, axis=0)
    gradmag = np.sqrt(gradrox**2 + gradroy**2 + gradroz**2)

    # Candidate cells: converging flow
    print('Identify candidate cells')
    temp = np.argwhere(divv < -convl)
    if ndim == 2:
        col, row, zrow = temp[:, 0], temp[:, 1], 0
    else:
        zrow, col, row = temp[:, 0], temp[:, 1], temp[:, 2]

    if not periodic:
        # Drop candidates whose +-avecyl sampling window would leave the domain
        w = max(int(avecyl), 2)
        ny_, nx_ = np.shape(divv)[-2:]
        keep = (col >= w) & (col <= ny_-1-w) & (row >= w) & (row <= nx_-1-w)
        if ndim == 3:
            nz_ = np.shape(divv)[0]
            keep = keep & (zrow >= w) & (zrow <= nz_-1-w)
            zrow = zrow[keep]
        print('Removing', int(np.sum(~keep)), 'candidates within', w, 'cells of the boundary')
        col = col[keep]
        row = row[keep]

    # One cell per front: keep only maxima of |grad rho| along the normal
    col, row, zrow = nonMaxSuppress(col, row, zrow, gradrox, gradroy, gradroz, gradmag,
                                    ndim, x, y, z, periodic)

    print('Begin loop for n=', np.size(col))
    time1 = time.perf_counter()
    args = (col, row, zrow, ds, gradrox, gradroy, gradroz, gradmag, divv, ndim, avecyl,
            shocktol, dbproject, nbrsearch, periodic)
    if nproc == 1:
        shocks = shockClassLoop(0, np.size(col)-1, *args)
    else:
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(nproc)
        print('using ', nproc, ' cores')
        sol = pool.starmap(shockClassLoop, [(a, b) + args for a, b in chunkRanges(np.size(col), nproc)])
        pool.close()
        pool.join()
        shocks = {k: np.vstack([s[k] for s in sol]) for k in sol[0]}
    time2 = time.perf_counter()
    print('Finding shocks took ', time2-time1, ' on ', nproc, ' cores')
    return shocks


###############################################################################
def chunkRanges(n, nproc):
    """Inclusive (start, end) index pairs splitting range(n) over nproc workers.
    Empty chunks come back as (0, -1) so the worker loop does nothing."""
    out = []
    for c in np.array_split(np.arange(n), nproc):
        out.append((int(c[0]), int(c[-1])) if c.size else (0, -1))
    return out


###############################################################################
def smoothdata(ro, vx, vy, vz, bx, by, bz, pr, smthfac):
    """Gaussian-smooth every variable with sigma = smthfac cells (0 leaves them unchanged).
    Scalar vz / bz (2D runs without out-of-plane components) are passed through.
    Returns the dict of variables used everywhere else: ro, vx, vy, vz, bx, by, bz, pr."""
    def sm(a):
        return gaussian_filter(a, smthfac) if np.size(a) > 1 else a
    return {'ro': sm(ro), 'vx': sm(vx), 'vy': sm(vy), 'vz': sm(vz), 'pr': sm(pr),
            'bx': sm(bx), 'by': sm(by), 'bz': sm(bz)}


###############################################################################
def divergence(ds, ndim, x, y, z):
    """div(v) from the coordinate arrays x, y (and z), so non-uniform grids are handled
    (np.gradient with coordinates: second order on stretched grids, the same stencil as a
    scalar spacing on uniform grids)."""
    if ndim == 2:
        divv = np.gradient(ds['vx'], x, axis=1)
        if np.size(ds['vy']) > 1:
            divv = divv + np.gradient(ds['vy'], y, axis=0)
        return divv
    return (np.gradient(ds['vx'], x, axis=2)
            + np.gradient(ds['vy'], y, axis=1)
            + np.gradient(ds['vz'], z, axis=0))


###############################################################################
def nonMaxSuppress(col, row, zrow, gradx, grady, gradz, gradmag, ndim, x, y, z, periodic=True):
    """Non-maximum suppression along the density-gradient normal (as in Canny edge detection).
    A candidate survives only if gradmag at the cell is >= gradmag one cell further along +n and
    > gradmag one cell back along -n (values bilinear/trilinear interpolated). The asymmetric
    tie-break keeps exactly one cell when the peak falls midway between two cells.
    The step is one cell in index space, so it also works on stretched grids.
    col, row, zrow are the y, x, z indices; returns the filtered (col, row, zrow)."""
    col = np.asarray(col, dtype=int)
    row = np.asarray(row, dtype=int)
    if np.size(col) == 0:
        return col, row, zrow
    mode = 'grid-wrap' if periodic else 'nearest'
    g0 = gradmag[(zrow, col, row) if ndim == 3 else (col, row)]
    ok = g0 > 1.0e-16
    safe = np.where(ok, g0, 1.0)
    dx = np.gradient(np.asarray(x, dtype=float))[row]
    dy = np.gradient(np.asarray(y, dtype=float))[col]
    if ndim == 2:
        ox = gradx[col, row]/safe/dx
        oy = grady[col, row]/safe/dy
        nrm = np.hypot(ox, oy)
        nrm[nrm == 0] = 1.0
        ox, oy = ox/nrm, oy/nrm
        gp = map_coordinates(gradmag, [col+oy, row+ox], order=1, mode=mode)
        gm = map_coordinates(gradmag, [col-oy, row-ox], order=1, mode=mode)
        keep = ok & (g0 >= gp) & (g0 > gm)
        print('Non-max suppression kept', int(keep.sum()), 'of', np.size(col), 'candidates')
        return col[keep], row[keep], zrow
    zrow = np.asarray(zrow, dtype=int)
    dz = np.gradient(np.asarray(z, dtype=float))[zrow]
    ox = gradx[zrow, col, row]/safe/dx
    oy = grady[zrow, col, row]/safe/dy
    oz = gradz[zrow, col, row]/safe/dz
    nrm = np.sqrt(ox**2 + oy**2 + oz**2)
    nrm[nrm == 0] = 1.0
    ox, oy, oz = ox/nrm, oy/nrm, oz/nrm
    gp = map_coordinates(gradmag, [zrow+oz, col+oy, row+ox], order=1, mode=mode)
    gm = map_coordinates(gradmag, [zrow-oz, col-oy, row-ox], order=1, mode=mode)
    keep = ok & (g0 >= gp) & (g0 > gm)
    print('Non-max suppression kept', int(keep.sum()), 'of', np.size(col), 'candidates')
    return col[keep], row[keep], zrow[keep]


###############################################################################
def shockClassLoop(istart, iend, col, row, zrow, ds, gradrox, gradroy, gradroz,
                   gradmag, divv, ndim, avecyl, tol=0.05, dbproject=False,
                   nbrsearch=0, periodic=True):
    """Classify candidate cells istart..iend (inclusive). Returns, per shock type k:
    k         : locations, [y,x] in 2D or [z,y,x] in 3D
    k+'mach'  : (slow, Alfven, fast) Mach numbers pre then post
    k+'qual'  : [mach_margin, bn_jump, bt_rotation_deg, shift], see shockQuality.
                shift = signed offset (cells along +n) of the cell whose classification
                was used; 0 = the candidate cell itself.
    k+'norm', k+'perp' : unit normal and tangent, (x, y[, z]) components, see shockVectors.
                Like mach and qual, they come from the cell whose classification was used.
    nbrsearch : if the candidate (a gradient maximum) cannot be classified, try the cells
                1..nbrsearch steps along its normal, nearest first (ties: larger Mach
                margin). The hit is still recorded at the candidate, so each front cell
                is counted at most once and fronts stay one cell thick."""
    hits = {k: [] for k in TRANS.values()}
    mach = {k: [] for k in TRANS.values()}
    qual = {k: [] for k in TRANS.values()}
    norm = {k: [] for k in TRANS.values()}
    perp = {k: [] for k in TRANS.values()}
    shape = np.shape(divv)
    w = max(int(avecyl), 2)

    def classifyCell(c, r, z):
        normarr = getNormVals(c, r, z, ds, gradrox, gradroy, gradroz, gradmag, ndim, avecyl)
        if dbproject:
            nnew = refineNormalDB(normarr, avecyl, ndim)
            if nnew is not None:
                normarr = getNormVals(c, r, z, ds, gradrox, gradroy, gradroz, gradmag, ndim,
                                      avecyl, normal=nnew)
        res = classifyShock(normarr, avecyl, tol=tol)
        if res is None:
            return None
        pre, pos, vsa, mpre, mpos = res
        key = TRANS.get((pre, pos))
        if key is None:
            return None
        return (key, list(mpre) + list(mpos), shockQuality(normarr, avecyl, mpre, mpos),
                shockVectors(normarr, avecyl, ndim))

    def neighbour(c, r, z, k):
        """Cell k steps along the gradient normal at (z,c,r), or None if outside the domain."""
        if ndim == 2:
            g = gradmag[c, r]
            if g <= 1e-16: return None
            idx = [c + k*gradroy[c, r]/g, r + k*gradrox[c, r]/g]
        else:
            g = gradmag[z, c, r]
            if g <= 1e-16: return None
            idx = [z + k*gradroz[z, c, r]/g, c + k*gradroy[z, c, r]/g, r + k*gradrox[z, c, r]/g]
        idx = [int(round(v)) for v in idx]
        dims = shape[-len(idx):]
        for j, (v, n) in enumerate(zip(idx, dims)):
            if periodic:
                idx[j] = v % n
            elif v < w or v > n-1-w:
                return None
        return (idx[0], idx[1], 0) if ndim == 2 else (idx[1], idx[2], idx[0])

    for i in range(int(istart), int(iend)+1):
        z = zrow if ndim == 2 else zrow[i]
        loc = [col[i], row[i]] if ndim == 2 else [zrow[i], col[i], row[i]]
        out = classifyCell(col[i], row[i], z)
        shift = 0
        k = 1
        while out is None and k <= nbrsearch:
            best = None
            for kk in (k, -k):
                nb = neighbour(col[i], row[i], z, kk)
                if nb is None: continue
                o = classifyCell(nb[0], nb[1], zrow if ndim == 2 else nb[2])
                if o is not None and (best is None or o[2][0] > best[0][2][0]):
                    best = (o, kk)
            if best is not None:
                out, shift = best
            k += 1
        if out is None:
            continue
        key, m, q, (nvec, tvec) = out
        hits[key].append(loc)
        mach[key].append(m)
        qual[key].append(q + [shift])
        norm[key].append(nvec)
        perp[key].append(tvec)
    shocks = {}
    for k in hits:
        shocks[k] = np.array(hits[k], dtype=int).reshape(-1, ndim)
        shocks[k+'mach'] = np.array(mach[k], dtype=float).reshape(-1, 6)
        shocks[k+'qual'] = np.array(qual[k], dtype=float).reshape(-1, 4)
        shocks[k+'norm'] = np.array(norm[k], dtype=float).reshape(-1, ndim)
        shocks[k+'perp'] = np.array(perp[k], dtype=float).reshape(-1, ndim)
    return shocks


###############################################################################
def getNormVals(col, row, zrow, var, gradx, grady, gradz, gradmag, ndim, avecyl,
                rad=False, normal=None):
    """Sample every variable along the shock normal through cell (zrow, col, row).

    The normal is the density-gradient direction at the cell, or `normal` (unit vector
    (nx, ny[, nz])) if given. 2*avecyl+1 points one cell apart are interpolated from the
    surrounding (2*avecyl+1)^ndim patch; the patch indices wrap periodically.
    rad=True also samples var['edref_m'] into 'rl'.

    Returns a dict of the sampled profiles (ro, vx, vy, [vz,] pr, bx, by, bz), the normal
    components normx, normy[, normz], and the normal / tangential components
    bperp, vperp (along n) and bpar, vpar (tangential), plus ang = angle of B to n.
    In 2D bpar/vpar are signed in-plane components (perpx, perpy is the in-plane tangent);
    in 3D they are the magnitudes of the tangential vectors."""
    if ndim == 2:
        g = gradmag[col, row]
        normx, normy = (0.0, 0.0) if g <= 1.0e-16 else (gradx[col, row]/g, grady[col, row]/g)
    else:
        g = gradmag[zrow, col, row]
        if g <= 1.0e-16:
            normx, normy, normz = 0.0, 0.0, 0.0
        else:
            normx, normy, normz = gradx[zrow, col, row]/g, grady[zrow, col, row]/g, gradz[zrow, col, row]/g
    if normal is not None:
        normx, normy = float(normal[0]), float(normal[1])
        if ndim == 3:
            normz = float(normal[2])

    # Patch indices around the cell, wrapped periodically
    shape = np.shape(var['ro'])
    off = np.arange(-avecyl, avecyl+1)
    tempx = (row + off) % shape[-1]
    tempy = (col + off) % shape[-2]

    nv = {}
    if ndim == 2:
        samp = lambda a: shocknormvals(a, tempx, tempy, normx, normy, avecyl)
        for k in ('ro', 'vx', 'vy', 'pr', 'bx', 'by'):
            nv[k] = samp(var[k])
        nv['bz'] = samp(var['bz']) if np.size(var['bz']) > 1 else var['bz']
        if rad:
            nv['rl'] = samp(var['edref_m'])
        nv['normx'] = normx
        nv['normy'] = normy
        sgn = 1.0 if normx >= 0 else -1.0
        nv['perpx'] = -normy*sgn
        nv['perpy'] = normx*sgn
        nv['bperp'] = nv['normx']*nv['bx'] + nv['normy']*nv['by']
        nv['bpar'] = nv['perpx']*nv['bx'] + nv['perpy']*nv['by']
        nv['vperp'] = nv['normx']*nv['vx'] + nv['normy']*nv['vy']
        nv['vpar'] = nv['perpx']*nv['vx'] + nv['perpy']*nv['vy']
        nv['ang'] = np.arctan(nv['bpar']/nv['bperp'])
        return nv

    tempz = (zrow + off) % shape[0]
    samp = lambda a: shocknormvals3d(a, tempx, tempy, tempz, normx, normy, normz, avecyl)
    for k in ('ro', 'vx', 'vy', 'vz', 'pr', 'bx', 'by'):
        nv[k] = samp(var[k])
    nv['bz'] = samp(var['bz']) if np.size(var['bz']) > 1 else var['bz']
    nv['normx'] = normx
    nv['normy'] = normy
    nv['normz'] = normz
    nv['bperp'] = nv['normx']*nv['bx'] + nv['normy']*nv['by'] + nv['normz']*nv['bz']
    nv['vperp'] = nv['normx']*nv['vx'] + nv['normy']*nv['vy'] + nv['normz']*nv['vz']
    b2 = nv['bx']**2 + nv['by']**2 + nv['bz']**2
    v2 = nv['vx']**2 + nv['vy']**2 + nv['vz']**2
    nv['bpar'] = np.sqrt(np.maximum(b2 - nv['bperp']**2, 0.0))   # |B_t|
    nv['vpar'] = np.sqrt(np.maximum(v2 - nv['vperp']**2, 0.0))   # |v_t|
    nv['ang'] = np.arctan2(nv['bpar'], np.abs(nv['bperp']))
    return nv


###############################################################################
def shocknormvals(var, tempx, tempy, normx, normy, avecyl):
    """Sample `var` (indexed [y,x]) at 2*avecyl+1 points one cell apart along the
    direction (normx, normy) through the centre of the patch tempy x tempx.

    tempx, tempy : integer x and y indices of the patch (already periodically wrapped)
    Bilinear interpolation inside the (2*avecyl+1)^2 patch with map_coordinates."""
    patch = var[np.ix_(np.asarray(tempy, dtype=int), np.asarray(tempx, dtype=int))]
    t = np.arange(-avecyl, avecyl+1, dtype=float)
    # line coordinates in patch index space; the centre cell is (avecyl, avecyl).
    # |n| <= 1, so every point lies inside the patch (mode only matters at round-off).
    coords = np.vstack([avecyl + normy*t, avecyl + normx*t])
    return map_coordinates(patch, coords, order=1, mode='nearest')


###############################################################################
def shocknormvals3d(var, tempx, tempy, tempz, normx, normy, normz, avecyl):
    """3D version of shocknormvals: `var` is indexed [z,y,x], the patch is
    tempz x tempy x tempx and the line follows (normx, normy, normz). Trilinear
    interpolation via map_coordinates."""
    patch = var[np.ix_(np.asarray(tempz, dtype=int), np.asarray(tempy, dtype=int),
                       np.asarray(tempx, dtype=int))]
    t = np.arange(-avecyl, avecyl+1, dtype=float)
    coords = np.vstack([avecyl + normz*t, avecyl + normy*t, avecyl + normx*t])
    return map_coordinates(patch, coords, order=1, mode='nearest')


###############################################################################
def refineNormalDB(normarr, avecyl, ndim, mindb=0.05):
    """Normal refined with div(B)=0: B_n is continuous, so the true normal is perpendicular
    to dB = B_post - B_pre. Returns the refined unit normal, or None if |dB| is too small
    (relative to |B| upstream) for the constraint to be meaningful."""
    ipre, ipos = prepostIndex(normarr['ro'], avecyl)
    comps = ('bx', 'by', 'bz') if ndim == 3 else ('bx', 'by')
    def val(k, i):
        v = np.asarray(normarr[k]); return float(v[i]) if v.size > 1 else float(v)
    dB = np.array([val(k, ipos) - val(k, ipre) for k in comps])
    Bmag = np.sqrt(sum(val(k, ipre)**2 for k in ('bx', 'by', 'bz')))
    nvec = np.array([normarr['normx'], normarr['normy']] + ([normarr['normz']] if ndim == 3 else []), dtype=float)
    if np.linalg.norm(dB) <= mindb*Bmag or np.linalg.norm(nvec) == 0:
        return None
    e = dB/np.linalg.norm(dB)
    nnew = nvec - (nvec @ e)*e
    mag = np.linalg.norm(nnew)
    if mag < 1e-3:              # dB parallel to n: constraint inconsistent, keep original
        return None
    return nnew/mag


###############################################################################
def classifyShock(normarr, avecyl, tol=0.05):
    """Classify one normal profile (output of getNormVals).
    Takes the pre/post states at the ends of the line, moves to the shock frame and computes
    the slow, Alfven and fast Mach numbers of the normal flow on each side.
    Returns (prestate, posstate, vsa, premach, posmach), or None if rejected (no compression,
    or the flow does not cross from the pre- to the post-shock side)."""
    ipre, ipos = prepostIndex(normarr['ro'], avecyl)
    vsa = getShockFrame(normarr['ro'][ipos], normarr['ro'][ipre],
                        normarr['vperp'][ipos], normarr['vperp'][ipre])
    if np.isnan(vsa):
        return None
    s = 1.0 if ipos > ipre else -1.0            # +1 if n points from pre to post
    if s*(normarr['vperp'][ipre] + vsa) <= 0:   # fluid not crossing pre -> post
        return None
    speeds = getWaveSpeeds(normarr['ro'], normarr['pr'], normarr['bx'],
                           normarr['by'], normarr['bz'], normarr['bperp'])
    states, machs = [], []
    with np.errstate(divide='ignore', invalid='ignore'):
        for i in (ipre, ipos):
            vf = np.abs(normarr['vperp'][i] + vsa)
            m = (vf/np.sqrt(speeds['vslow2'][i]),
                 vf/np.sqrt(speeds['vap2'][i]),
                 vf/np.sqrt(speeds['vfast2'][i]))
            machs.append(m)
            states.append(getState(*m, tol=tol))
    return states[0], states[1], vsa, machs[0], machs[1]


###############################################################################
def shockQuality(normarr, avecyl, mpre, mpos):
    """[mach_margin, bn_jump, bt_rotation_deg] for a classified shock.
    mach_margin : smallest |M-1| over the six finite Mach numbers (how close the
                  classification is to flipping)
    bn_jump     : |Bn_post - Bn_pre| / |B_pre|; should be ~0 (div B = 0), large values
                  mean the normal is unreliable. With dbproject=True the refined normal
                  makes this ~0 by construction, so it is only informative without it
                  (or where the refinement was skipped because |dB| was small).
    bt_rotation : angle (deg) between the tangential field before and after the front.
                  Coplanar slow/fast shocks: ~0; 2-4 / 1-3 intermediate shocks: ~180."""
    m = np.abs(np.r_[mpre, mpos] - 1.0)
    m = m[np.isfinite(m)]
    margin = float(np.min(m)) if m.size else np.nan
    ipre, ipos = prepostIndex(normarr['ro'], avecyl)
    bz = np.asarray(normarr['bz'])
    bzpre = float(bz[ipre]) if bz.size > 1 else float(bz)
    Bmag = np.sqrt(normarr['bx'][ipre]**2 + normarr['by'][ipre]**2 + bzpre**2)
    bnj = abs(normarr['bperp'][ipos] - normarr['bperp'][ipre])/Bmag if Bmag > 1e-30 else np.nan
    n = np.array([normarr['normx'], normarr['normy'], normarr.get('normz', 0.0)], dtype=float)
    def bt(i):
        bzi = float(bz[i]) if bz.size > 1 else float(bz)
        B = np.array([normarr['bx'][i], normarr['by'][i], bzi])
        return B - (B @ n)*n
    a, b = bt(ipre), bt(ipos)
    den = np.linalg.norm(a)*np.linalg.norm(b)
    rot = float(np.degrees(np.arccos(np.clip(a @ b/den, -1, 1)))) if den > 1e-30 else np.nan
    return [margin, float(bnj), rot]


###############################################################################
def shockVectors(normarr, avecyl, ndim, mindb=1.0e-3, mintan=0.05):
    """Unit normal and tangent of a classified front, as (x, y[, z]) components.

    norm : the normal used to sample the profile (after any dbproject refinement),
           flipped if needed so it points from the pre-shock to the post-shock sample.
    perp : 2D - the in-plane tangent z_hat x norm = (-ny, nx), so (norm, perp) is
                right-handed. (The internal perpx, perpy used for bpar/vpar flips sign
                with normx instead.)
           3D - direction of the jump in tangential field, dB_t = dB - (dB.n)n with
                dB = B_post - B_pre. By div(B)=0 the jump lies in the front, so this is
                robust to normal errors, and for a coplanar shock it spans the plane of
                B and n. Sign: the upstream B_t lies along +perp for fast shocks (B_t
                grows) and along -perp for slow and intermediate shocks.
                Only if |dB_t| <= mindb*|B_pre| (no tangential jump: hydro-like or parallel
                shocks) does it fall back to the upstream B_t (if > mintan*|B_pre|), then
                the upstream tangential velocity, then any fixed unit vector perpendicular
                to n; the sign rule above then does not apply.
                The second tangent is np.cross(norm, perp)."""
    ipre, ipos = prepostIndex(normarr['ro'], avecyl)
    comps = ('normx', 'normy', 'normz')[:ndim]
    n = np.array([normarr[k] for k in comps], dtype=float)
    mag = np.linalg.norm(n)
    if mag > 0:
        n = n/mag
    if ipos < ipre:                       # make n point from pre- to post-shock side
        n = -n
    if ndim == 2:
        return n.tolist(), [-n[1], n[0]]

    def vec(prefix, i):
        out = []
        for c in 'xyz':
            v = np.asarray(normarr[prefix + c])
            out.append(float(v[i]) if v.size > 1 else float(v))
        return np.array(out)

    Bpre, Bpos, vpre = vec('b', ipre), vec('b', ipos), vec('v', ipre)
    Bref = np.linalg.norm(Bpre)
    for cand, ref, tol in ((Bpos - Bpre, Bref, mindb), (Bpre, Bref, mintan),
                           (vpre, np.linalg.norm(vpre), mintan)):
        tv = cand - (cand @ n)*n
        tm = np.linalg.norm(tv)
        if tm > tol*ref and tm > 1.0e-30:
            return n.tolist(), (tv/tm).tolist()
    # nothing tangential: any unit vector perpendicular to n
    axis = np.eye(3)[int(np.argmin(np.abs(n)))]
    tv = axis - (axis @ n)*n
    return n.tolist(), (tv/np.linalg.norm(tv)).tolist()


###############################################################################
def prepostIndex(ro, avecyl, off=None):
    """Indices of the pre- and post-shock samples along the normal profile.

    ro     : density sampled along the normal, length 2*avecyl+1, centre = avecyl
    off    : distance in samples from the centre (default avecyl, i.e. the window ends)
    Returns (ipre, ipos); the post-shock side is the denser one."""
    off = avecyl if off is None else off
    lo, hi = avecyl - off, avecyl + off       # symmetric about the front
    if ro[lo] < ro[hi]:
        return lo, hi                         # (ipre, ipos): density rises along n
    return hi, lo                             # profile reversed: swap


###############################################################################
def getShockFrame(ropos, ropre, vperppos, vperppre, tol=1.0e-2):
    """Shock-frame offset vsa = -Vs from mass conservation, rho_pre*(v_pre+vsa) =
    rho_post*(v_post+vsa). Returns NaN unless the post/pre density ratio exceeds 1+tol
    (no abs: expansions fail)."""
    drho = ropos - ropre
    if (not np.isfinite(drho)) or drho <= tol*ropre:
        return np.nan
    return -(ropos*vperppos - ropre*vperppre)/drho


###############################################################################
def getWaveSpeeds(ro, pr, bx, by, bz, bperp, gamma=5.0/3.0):
    """Phase speeds squared of the slow/Alfven/fast waves along the shock normal.

    bperp : component of B along the shock normal (signed or unsigned)
    bz    : may be a scalar 0 in 2D runs without a z-component
    Returns a dict with cs2, va2 (full |B|), vap2 (along n), vslow2, vfast2."""
    b2   = bx**2 + by**2 + bz**2                 # full |B|^2, including bz
    cs2  = gamma*pr/ro
    va2  = b2/ro
    vap2 = bperp**2/ro                           # Alfven speed along n

    # cos^2(theta_Bn) = Bn^2/|B|^2 ; no tangential direction needed
    b2 = np.broadcast_to(b2, np.shape(ro))
    cos2 = np.divide(np.broadcast_to(bperp, np.shape(ro))**2, b2,
                     out=np.ones(np.shape(ro)), where=b2 > 1e-30)

    disc = np.sqrt(np.maximum((cs2+va2)**2 - 4.0*va2*cs2*cos2, 0.0))
    return {'cs2': cs2, 'va2': va2, 'vap2': vap2,
            'vslow2': 0.5*(cs2+va2-disc),
            'vfast2': 0.5*(cs2+va2+disc)}


###############################################################################
def getState(vslow, valp, vfast, tol=0.0):
    """State of the flow from its slow, Alfven and fast Mach numbers:
    1: M_fast>=1 | 2: M_alf>=1>M_fast | 3: M_slow>=1>M_alf | 4: all <1
    0: marginal (any Mach number within tol of 1) or undefined (NaN)."""
    m = np.array([vslow, valp, vfast], dtype=float)
    if np.any(np.isnan(m)) or np.any(np.abs(m - 1.0) < tol):
        return 0
    if vfast >= 1: return 1
    if valp  >= 1: return 2
    if vslow >= 1: return 3
    return 4


###############################################################################
def restoreShocks(fname):
    """Load a shocks dict previously saved to HDF5 (one dataset per key)."""
    import h5py
    with h5py.File(fname, "r") as hf:
        return {f: np.asarray(hf[f]) for f in hf.keys()}


###############################################################################
def shockFilter(shocks, maxDis):
    """Remove isolated shock points: a point is dropped if its nearest other point of the
    same type is at least maxDis cells away. Works for 2D [y,x] and 3D [z,y,x] locations.
    Every shock type is returned (possibly unchanged), with matching '<type>mach',
    '<type>qual', '<type>norm' and '<type>perp' rows."""
    from scipy.spatial import cKDTree
    shocks2 = {}
    for k in TRANS.values():
        if k not in shocks:
            continue
        pts = np.asarray(shocks[k])
        keep = np.ones(len(pts), dtype=bool)
        if len(pts) >= 3:
            d, _ = cKDTree(pts).query(pts, k=2)    # d[:,1] = distance to nearest other point
            keep = d[:, 1] < maxDis
        shocks2[k] = pts[keep].astype('int')
        for suffix in ('mach', 'qual', 'norm', 'perp'):
            if k+suffix in shocks:
                shocks2[k+suffix] = np.asarray(shocks[k+suffix])[keep]
    return shocks2


###############################################################################
def shockLine(loc, ds, avecyl=5, ndim=2, rad=False, getPoints=False, getEnergy=False):
    """Diagnostics for one shock location (e.g. a row of shocks['slow']).

    loc : [y, x] (2D) or [z, y, x] (3D) indices, as returned by shockid()
    ds  : dict with ro, vx, vy, vz, bx, by, bz, pr and coordinate arrays xgrid, ygrid
          (and zgrid in 3D). Uses the data as given (no smoothing) and the raw density
          gradient for the normal.
    Returns a dict with the normal profiles 'normarr' (see getNormVals), 'ipre', 'ipos',
    the shock-frame offset 'vsa' and the wave 'speeds' along the line. Optional:
      getPoints : 'normlinex', 'normliney', the x and y indices of the sampled points
      getEnergy : kinetic 'ke', magnetic 'me' and thermal 'te' energy density profiles"""
    if ndim == 3:
        zrow, col, row = loc[0], loc[1], loc[2]
        z = ds['zgrid']
    else:
        zrow, col, row = 0, loc[0], loc[1]
        z = 0.0
    x = ds['xgrid']
    y = ds['ygrid']

    gradroz = 0.0
    if ndim == 2:
        gradrox = np.gradient(ds['ro'], x, axis=1)
        gradroy = np.gradient(ds['ro'], y, axis=0)
    else:
        gradrox = np.gradient(ds['ro'], x, axis=2)
        gradroy = np.gradient(ds['ro'], y, axis=1)
        gradroz = np.gradient(ds['ro'], z, axis=0)
    gradmag = np.sqrt(gradrox**2 + gradroy**2 + gradroz**2)

    normarr = getNormVals(col, row, zrow, ds, gradrox, gradroy, gradroz, gradmag, ndim, avecyl,
                          rad=rad)
    ipre, ipos = prepostIndex(normarr['ro'], avecyl)
    vsa = getShockFrame(normarr['ro'][ipos], normarr['ro'][ipre],
                        normarr['vperp'][ipos], normarr['vperp'][ipre])
    speeds = getWaveSpeeds(normarr['ro'], normarr['pr'], normarr['bx'], normarr['by'],
                           normarr['bz'], normarr['bperp'])
    lineData = {'vsa': vsa, 'speeds': speeds, 'normarr': normarr, 'ipre': ipre, 'ipos': ipos}

    if getPoints:
        t = np.linspace(-avecyl, avecyl, 2*avecyl+1)
        lineData['normlinex'] = normarr['normx']*t + row     # x index
        lineData['normliney'] = normarr['normy']*t + col     # y index

    if getEnergy:
        lineData['ke'] = 0.5*normarr['ro']*(normarr['vpar']**2 + normarr['vperp']**2)
        lineData['me'] = 0.5*(normarr['bpar']**2 + normarr['bperp']**2)
        lineData['te'] = normarr['pr']/(5.0/3.0 - 1.0)
    return lineData


###############################################################################
def xyzTozyx(arr):
    """Reorder a 3D array from [x,y,z] to the [z,y,x] layout shockid() expects."""
    return np.transpose(arr, (2, 1, 0))
