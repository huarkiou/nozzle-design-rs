import scipy.optimize as so
import scipy.integrate as si
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as pts


def cal_cutellipse_area(a: float, b: float, x0: float) -> float:
    """ 计算半椭圆截短后部分的面积 x0为沿长半轴a的截短位置到长半轴中间对称轴的距离 面积为此对称轴到截短位置之间的椭圆面积 """
    if x0 > a:
        x0 = a
        print("Warning: x0 should lower than a")
    elif x0 < 0:
        x0 = 0
        print("Warning: x0 should greater than 0")
    return si.quad(lambda x: 2 * b * math.sqrt(1 - (x / a)**2), 0, x0)[0]


def cal_halfellipse_area(ellipse: tuple[float, float]) -> float:
    return cal_cutellipse_area(ellipse[0], ellipse[1], ellipse[0])


def main() -> int:
    write_points_by_params(
        (0.9, 0.6), (30, 0.6), (-0.05, 0.1),
        r"D:\Projects\CFD\nozzledual2\model01\shapepoints.txt")
    write_points_by_params(
        (0.9, 0.6), (1.5, 0.6), (-0.05, 0.1),
        r"D:\Projects\CFD\nozzledual2\model02\shapepoints.txt")
    write_points_by_params(
        (0.9, 0.6), (0.905, 0.6), (-0.05, 0.1),
        r"D:\Projects\CFD\nozzledual2\model03\shapepoints.txt")

    return 0


def write_points_by_params(base: tuple[float, float], greater: tuple[float,
                                                                     float],
                           offset: tuple[float, float], filepath: str) -> None:
    base_ellipse = base
    greater_ellipse = greater
    x0, y0 = offset

    area0 = cal_halfellipse_area(base_ellipse)
    xt = so.fsolve(
        lambda x: area0 - cal_cutellipse_area(greater_ellipse[
            0], greater_ellipse[1], x), greater_ellipse[0] / 2)[0]
    theta = math.acos(xt / greater_ellipse[0])
    print("Setting points in file ", filepath)
    print("Input: ", base, greater, offset)
    print(f"xt={xt} θ=theta({math.degrees(theta)})°")

    # 组装三者 1.切割后的半椭圆 2.直线 3.半椭圆
    x = np.arange(0, 1, 0.02) * xt + x0
    x = np.append(x, xt + x0)
    y = np.zeros(x.shape)
    for i in range(0, y.shape[0]):
        y[i] = y0 + greater_ellipse[1] * math.sqrt(1 - (
            (x[i] - x0) / greater_ellipse[0])**2)
    x_tmp = np.flip(x)
    y_tmp = np.zeros(x_tmp.shape)
    for i in range(0, y_tmp.shape[0]):
        y_tmp[i] = y0 - greater_ellipse[1] * math.sqrt(1 - (
            (x_tmp[i] - x0) / greater_ellipse[0])**2)
    x = np.append(x, x_tmp)
    y = np.append(y, y_tmp)

    t_tmp = np.arange(math.pi / 2, math.pi / 2 * 3, 0.1)
    x_tmp = np.zeros(t_tmp.shape)
    y_tmp = np.zeros(t_tmp.shape)
    for i in range(0, t_tmp.shape[0]):
        x_tmp[i] = x0 + base_ellipse[0] * math.cos(t_tmp[i])
        y_tmp[i] = y0 + base_ellipse[1] * math.sin(t_tmp[i])
    x = np.append(x, np.flip(x_tmp))
    y = np.append(y, np.flip(y_tmp))

    x[np.abs(x) < np.finfo(np.float32).eps] = 0
    y[np.abs(y) < np.finfo(np.float32).eps] = 0

    result = np.vstack((x, y)).T
    np.savetxt(filepath, result, delimiter=' ')

    ## 绘图
    fig = plt.figure()
    # 参考单位圆
    refcircle = pts.Circle(xy=(0, 0), radius=1., fill=False, linestyle='--')
    plt.gca().add_patch(refcircle)
    refellipse = pts.Ellipse(xy=offset,
                             width=2 * greater[0],
                             height=2 * greater[1],
                             fill=False)
    plt.gca().add_patch(refellipse)

    plt.axis("equal")
    plt.plot(x, y)
    plt.show()


if __name__ == "__main__":
    ret = main()
    exit(ret)
