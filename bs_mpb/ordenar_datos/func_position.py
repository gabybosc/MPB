import numpy as np


def boundary_position(t_boundary, t, x, y, z):
    """
    Finds spacecraft cartesian position vector at the time of a boundary
    (Bs or MPB) crossing.

    Parameters:

        t_boundary: float
                    crossing time of boundary in decimal hour format

        t: N-1 array
          time vector in decimal hour format

        x,y,z: N-1 arrays
              spacecraft position vector components


      Returns:

          r_boundary: N-3 array
                     spacecraft position vector at crossing time
    """

    i_boundary = (abs(t_boundary - t)).argmin()

    x_boundary = x[i_boundary]
    y_boundary = y[i_boundary]
    z_boundary = z[i_boundary]

    r_boundary = np.array([x_boundary, y_boundary, z_boundary]).T

    return r_boundary


def cartesian2polar(x, y, z, r0, sym_axis="x"):
    """
    Obtain polar coordinates from cartesian components.

    Parameters:
        x,y,z: float
               position of 1 point in conic surface measured from origin

        r0: (1,3) array
            focus of conic surface

        sym_axis: str
                  axis of symmetry for conic surface (e.g.: 'x')

    Returns:
        r:  float
            distance of point in conic surface measured from focus

        theta: float
               polar angle measured from symmetry axis

        phi: float
             azimuth angle measured in plane perpendicular to symmetry axis
    """

    X = x - r0[0]
    Y = y - r0[1]
    Z = z - r0[2]

    r = np.sqrt(X**2 + Y**2 + Z**2)

    if sym_axis == "x":
        theta = np.arccos(X / r)
        phi = np.arctan2(Z, Y)

    elif sym_axis == "y":
        theta = np.arccos(Y / r)
        phi = np.arctan2(X, Z)

    elif sym_axis == "z":
        theta = np.arccos(Z / r)
        phi = np.arctan2(Y, X)

    return r, theta, phi  # type: ignore
