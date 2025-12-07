import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scienceplots

plt.style.use(["science", "nature"])

plt.rcParams.update(
    {
        "font.family": "serif",  # specify font family here
        "font.serif": ["Times New Roman"],  # specify font here
        # "font.size": 11
    }
)  # specify font size here

# mpl.use("pgf")  # 修改绘图后端
pgf_with_rc_fonts = {
    "font.family": ["Times New Roman", "SimSun"],
    "font.serif": ["Times New Roman", "SimSun"],
    "font.sans-serif": ["Times New Roman", "SimSun"],
    "mathtext.fontset": "custom",
}
mpl.rcParams.update(pgf_with_rc_fonts)


def main():
    for a in range(-9, 10, 3):
        plot_line(a)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.legend()
    plt.show()
    # plt.savefig(r"fig1.png", dpi=300)


def plot_line(a: float):
    x = np.arange(0, 1.01, 0.01)
    y = x
    if a > 0:
        y = np.arctan((2 * x - 1) * a) / np.arctan(a) / 2.0 + 0.5
    elif a < 0:
        y = (np.tan(np.arctan(a) * (2 * x - 1)) / a + 1) / 2.0
    else:
        y = x
    plt.plot(x, y, label=r"$a={:+2}$".format(a))


if __name__ == "__main__":
    main()
