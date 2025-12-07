"""
    从指定文件中读取空间坐标点，拟合椭圆方程
"""

import matplotlib.pyplot as plt
import scipy.optimize as so
import numpy as np


def main():
    sourcefile = r"D:\Z_Program\nozzleproject\t240705\点集\inlet_parallel.csv"
    data = np.loadtxt(sourcefile, delimiter=',')

    x0 = data[:, 0].mean()
    y0, z0, r, res_optimized = fit_ellipse(data[:, 1],
                                           data[:, 2],
                                           r'$y$',
                                           r'$z$',
                                           tol=1e-13)
    print(f"message: {res_optimized.message}\njac: {res_optimized.jac}\n")

    print(f"圆心坐标:({x0}, {y0}, {z0}) \n半径r = {r}")


def my_fun(parameters, x_samples, y_samples):
    # 两焦点坐标以及圆上的点到两焦点的距离的和作为优化参数
    x0, y0, sum_of_target_distance_between_edge_and_center = parameters
    # 计算实际距离
    sum_of_actual_distance_between_edge_and_center = ((x_samples - x0)**2 +
                                                      (y_samples - y0)**2)**0.5

    # print(np.average(sum_of_actual_distance_between_edge_and_two_focus))
    # 返回方差
    return np.sum(((sum_of_actual_distance_between_edge_and_center -
                    sum_of_target_distance_between_edge_and_center)**2) /
                  (len(x_samples) - 1))


def fit_ellipse(x_samples, y_samples, xlabel, ylabel, tol=1e-6):
    # 归一化
    vmax = max(np.max(x_samples), np.max(y_samples))
    # x_samples = x_samples / vmax
    # y_samples = y_samples / vmax
    # 优化
    res_optimized = so.minimize(fun=my_fun,
                                x0=np.array([0., 0., 1.]),
                                args=(x_samples, y_samples),
                                tol=tol)

    x0_res, y0_res, l2_res = res_optimized.x
    # 依据优化得到的函数生成椭圆曲线
    # 计算半径
    r_res = l2_res
    # 计算椭圆中心坐标
    x0 = x0_res
    y0 = y0_res

    # 极坐标轴序列
    theta_res = np.linspace(0.0, 6.28, 1000)
    # 生成椭圆上的点
    x_res = x0 + r_res * np.cos(theta_res)
    y_res = y0 + r_res * np.sin(theta_res)

    # plt.style.use("one")
    plt.axis('equal')
    plt.scatter(x_samples,
                y_samples,
                color="magenta",
                marker="+",
                zorder=1,
                s=80,
                label="samples")
    plt.plot(x_res, y_res, color="deepskyblue", zorder=2, label="fitted curve")
    plt.scatter(np.array([x0_res]),
                np.array([y0_res]),
                zorder=3,
                color="r",
                label="focus point")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    # plt.savefig("Figsave/a={:.3f};b={:.3f};theta={:.2f}deg.svg".format(a_res, b_res, alpha_res))
    plt.show()
    return x0, y0, r_res, res_optimized


if __name__ == "__main__":
    main()
