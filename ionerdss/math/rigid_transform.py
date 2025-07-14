import numpy as np

def rigid_transform_3d(P, Q):
    """
    Compute the optimal rigid body transform (rotation R and translation t)
    that aligns two 3D point sets P and Q.

    Parameters
    ----------
    P : ndarray of shape (N, 3)
        Source point cloud.
    Q : ndarray of shape (N, 3)
        Target point cloud.

    Returns
    -------
    R : ndarray of shape (3, 3)
        Rotation matrix.
    t : ndarray of shape (3,)
        Translation vector.
    """
    P = np.asarray(P)
    Q = np.asarray(Q)

    assert P.shape == Q.shape, "Point sets must be the same shape"

    P_mean = P.mean(axis=0)
    Q_mean = Q.mean(axis=0)

    P_centered = P - P_mean
    Q_centered = Q - Q_mean

    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Correct for reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = Q_mean - R @ P_mean
    return R, t


def apply_transform(R, t, X):
    """
    Apply a rigid transformation (rotation + translation) to a point or set of points.

    Parameters
    ----------
    R : ndarray of shape (3, 3)
        Rotation matrix.
    t : ndarray of shape (3,)
        Translation vector.
    X : ndarray of shape (3,) or (N, 3)
        Point(s) to transform.

    Returns
    -------
    ndarray of shape like X
        Transformed coordinates.
    """
    X = np.asarray(X)
    return (R @ X.T).T + t