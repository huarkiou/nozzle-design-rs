import matplotlib.pyplot as plt
from math import cos, sin, sqrt

p1: tuple[float, float, float] = (1.8951087644550260, 0.0000000000000000, 0.0)
p2: tuple[float, float, float] = (
    1.8919379416868118,
    0.0013837169816496176,
    0.00052976493774764544,
)
pr: tuple[float, float, float] = (1.898279448151052, 0.0013840356254392625, 0)


MAX_DISTANCE: float = 1.0


def distance(p1: tuple[float, float, float], p2: tuple[float, float, float]) -> float:
    return sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


MAX_DISTANCE = max(distance(p1, p2), distance(p1, pr), distance(p2, pr))


def main():
    plt.figure()
    plot_point_curve(p1)
    plot_point_curve(p2)
    plot_point_curve(pr)
    plt.show()


def plot_point_curve(point: tuple[float, float, float]):
    plt.plot(point[0], point[1], "o")
    theta = point[2]
    length: float = MAX_DISTANCE / 2.0
    p1 = (point[0] - length * cos(theta), point[1] - length * sin(theta))
    p2 = (point[0] + length * cos(theta), point[1] + length * sin(theta))
    plt.plot((p1[0], p2[0]), (p1[1], p2[1]), "-")


if __name__ == "__main__":
    main()
