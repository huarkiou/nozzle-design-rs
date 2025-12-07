import matplotlib.pyplot as plt
import matplotlib.patches as pts
import numpy as np


def main():
    fig = plt.figure()

    # 参考单位圆
    refcircle = pts.Circle(xy=(0, 0), radius=1., fill=False, linestyle='--')
    plt.gca().add_patch(refcircle)

    # 进口截面
    secin = pts.Circle(xy=(0.0, 0.4), radius=0.5, fill=False, color='red')
    plt.gca().add_patch(secin)

    # 出口截面
    secout = pts.Rectangle(xy=(0.0, 0.35),
                           width=1.2,
                           height=0.8,
                           angle=0,
                           rotation_point='center',
                           fill=False,
                           color='blue')
    secout.set_xy((secout.get_x() - secout.get_width() / 2.,
                   secout.get_y() - secout.get_height() / 2.))
    plt.gca().add_patch(secout)

    # 显示
    plt.axis("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


if __name__ == "__main__":
    main()
