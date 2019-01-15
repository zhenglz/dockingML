import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math

def fitPlane(points) :
    """
    fit some points to a plane
    https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
    plane function: ax + by - z + c = 0
    or: ax + by - z = -c
    Ax + By + Cz + D = 0
    we need to determine [a, b, c]
    A B C D = a  b -1 c
    :param points: ndarray, M*3
    :return: array, [a, b, c]
    """

    # prepare dataset
    xs = np.array(points)[:, 0]
    ys = np.array(points)[:, 1]
    zs = np.array(points)[:, 2]

    tmp_A = []
    tmp_B = []
    for i in range(xs.shape[0]):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_B.append(zs[i])

    B = np.matrix(tmp_B).T
    A = np.matrix(tmp_A)

    # do fit
    fit = (A.T * A).I * A.T * B
    errors = B - A * fit
    residual = np.linalg.norm(errors)

    return fit

def point_distance(params, point) :
    """
    determine the distance between a point to a plane
    :param params:
    :param point:
    :return:
    """

    distance = math.sqrt(params[0] ** 2 + params[1] ** 2 + 1)

    distance = (params[0] * point[0] +
                params[1] * point[1] +
                (-1.0) * point[2] +
                params[2]
                ) / distance

    return distance

if __name__ == "__main__" :

    points = [
        [0.0, 0.0, 0.1 ],
        [1.0, 1.0, 0.0 ],
        [1.0, 0.0, -0.2 ],
        [0.5, 0.5, 0.0 ],
        [-1.0, 2.0, 0],
        [0.5, -0.5, 0.2],
        [0.2, 0.2, 0.01],
        [0.7, 0.9, 0.002],
        [4.5, 2.3, -0.1],
        [-2.5, -5.0, 0.02]
    ]

    p = [0.0, 1.0, 5.0 ]

    fit = fitPlane(points)

    points.append(p)

    xs = np.array(points)[:, 0]
    ys = np.array(points)[:, 1]
    zs = np.array(points)[:, 2]

    d = point_distance(fit, p)
    print("Distance is %f" % d)

    # plot raw data
    plt.figure()
    ax = plt.subplot(111, projection='3d')
    ax.scatter(xs, ys, zs, color='orange')

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    X, Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                       np.arange(ylim[0], ylim[1]))
    Z = np.zeros(X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r, c] = fit[0] * X[r, c] + fit[1] * Y[r, c] + fit[2]
    ax.plot_wireframe(X, Y, Z, color='k')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

