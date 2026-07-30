#It reads phase velocities from two overlapping subarray inversions, 
#then replaces duplicate grid nodes in the overlap zone with a Gaussian distance-weighted average
 #giving more trust to whichever subarray's center is closer
awk '
{
    key = $1 FS $2

    # Save the order of first appearance
    if (!(key in seen)) {
        seen[key] = 1
        order[++m] = key
        xlon[key] = $1
        xlat[key] = $2
    }

    v[key, ++n[key]] = $3
}
END{
    cW = -118.5    # center of western array
    cE = -109.5    # center of eastern array
    L  = 4.0       # characteristic length (degrees)

    # Loop in original input order
    for (i = 1; i <= m; i++) {
        key = order[i]
        xlo = xlon[key] + 0
        xla = xlat[key] + 0

        if (n[key] == 1) {
            val = v[key,1]
        }
        else {
            vW = v[key,1]
            vE = v[key,2]

            dW = xlo - cW
            if (dW < 0) dW = -dW

            dE = xlo - cE
            if (dE < 0) dE = -dE

            wW = exp(-(dW/L)^2)
            wE = exp(-(dE/L)^2)

            val = (wW*vW + wE*vE) / (wW + wE)
        }

        printf "%12.6f %12.7f %14.8f\n", xlo, xla, val
    }
}
' overlapeast_west.d > smooth_overlapeastwest.d