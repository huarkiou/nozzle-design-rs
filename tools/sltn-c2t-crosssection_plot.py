import matplotlib.pyplot as plt
import matplotlib.patches as pts
import numpy as np


def main():
    fig = plt.figure()

    # 参考单位圆
    refcircle = pts.Circle(xy=(0, 0), radius=1., fill=False, linestyle='--')
    plt.gca().add_patch(refcircle)

    # 进口截面
    secin = pts.Circle(xy=(0.0, 0.25), radius=0.75, fill=False, color='red')
    plt.gca().add_patch(secin)

    # 出口截面
    secout = pts.RegularPolygon(xy=(0.0, 0.25),
                                numVertices=3,
                                radius=0.73,
                                fill=False,
                                color='blue')
    # print(secout.get_path().vertices)
    plt.gca().add_patch(secout)
    # 显示
    plt.axis("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
