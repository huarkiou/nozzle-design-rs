import numpy as np
import scipy.optimize as so


def main():
    sourcefile = r"D:\Z_Program\nozzleproject\t240705\点集\inlet.csv"
    data = np.loadtxt(sourcefile, delimiter=',')

    p = estimate_plane_with_leastsq(data)
    center, normal, origin, plane_x_pt, plane_y_pt = get_proper_plane_params(
        p, data)
    print("estimate plane info:"
          f"  center: {center}\n"
          f"  normal: {normal}\n"
          f"  origin: {origin}\n"
          f"  plane_x_pt: {plane_x_pt}\n"
          f"  plane_y_pt: {plane_y_pt}")


def func_plane(p, x, y):
    """ 数据拟合函数 """
    a, b, c = p
    return a * x + b * y + c


def residuals(p, x, y, z):
    """ 误差函数 """
    return z - func_plane(p, x, y)


def estimate_plane_with_leastsq(pts: np.ndarray):
    """ 根据最小二乘拟合出平面参数 """
    p0 = [1, 0, 1]
    plsq = so.leastsq(residuals, p0, args=(pts[:, 0], pts[:, 1], pts[:, 2]))
    return plsq[0]


def get_proper_plane_params(p, pts: np.ndarray):
    """ 根据拟合的平面的参数，得到实际显示的最佳的平面 """
    np_pts_mean = np.mean(pts, axis=0)
    np_pts_min = np.min(pts, axis=0)
    np_pts_max = np.max(pts, axis=0)

    plane_center_z = func_plane(p, np_pts_mean[0], np_pts_mean[1])
    plane_center = [np_pts_mean[0], np_pts_mean[1], plane_center_z]

    plane_origin_z = func_plane(p, np_pts_min[0], np_pts_min[1])
    plane_origin = [np_pts_min[0], np_pts_min[1], plane_origin_z]

    if np.linalg.norm(p) < 1e-10:
        print(r'plsq 的 norm 值为 0 {}'.format(p))
    plane_normal = p / np.linalg.norm(p)

    plane_pt_1 = [
        np_pts_max[0], np_pts_min[1],
        func_plane(p, np_pts_max[0], np_pts_min[1])
    ]
    plane_pt_2 = [
        np_pts_min[0], np_pts_max[1],
        func_plane(p, np_pts_min[0], np_pts_max[1])
    ]
    return plane_center, plane_normal, plane_origin, plane_pt_1, plane_pt_2


if __name__ == '__main__':
    main()
