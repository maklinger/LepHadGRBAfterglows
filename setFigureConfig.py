import matplotlib.pyplot as plt
import numpy as np
import os

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

def configure(config):
    alpha = config["alpha"]
    plt.rc('font', size=alpha*config["small"])
    plt.rc('axes', titlesize=alpha*config["small"])
    plt.rc('axes', labelsize=alpha*config["medium"])
    plt.rc('xtick', labelsize=alpha*config["small"])
    plt.rc('ytick', labelsize=alpha*config["small"])
    plt.rc('legend',
           fontsize=alpha*config["small"],
           handlelength=config["legendHandleLength"]
           )
    plt.rc('figure',
           titlesize=alpha*config["big"],
           figsize=(config["width"], config["aspRatio"]*config["width"])
           )


# a plot with DinA4 full page width
def fullWidth(alpha=1, aspectRatio=4./6., legendHandleLength=2):
    config = {}
    config["alpha"] = alpha
    config["aspRatio"] = aspectRatio

    config["small"] = 12
    config["medium"] = 14
    config["big"] = 16

    config["width"] = 6.5

    config["legendHandleLength"] = legendHandleLength

    configure(config)


# a plot with DinA4 half page width
def halfWidth(alpha=1,  aspectRatio=4./6., legendHandleLength=1):
    config = {}
    config["alpha"] = alpha
    config["aspRatio"] = aspectRatio

    config["small"] = 7
    config["medium"] = 8
    config["big"] = 9

    config["width"] = 3.25

    config["legendHandleLength"] = legendHandleLength

    configure(config)

# a plot for Jupyter notebook
def jupyter(alpha=1.3,  aspectRatio=4./6., legendHandleLength=1):
    config = {}
    config["alpha"] = alpha
    config["aspRatio"] = aspectRatio

    config["small"] = 7
    config["medium"] = 8
    config["big"] = 9

    config["width"] = 10

    config["legendHandleLength"] = legendHandleLength

    configure(config)


# save figure as pdf and png in one command
def saveFigure(fig, name, pdf=True, png=True, dpi=400):
    if pdf:
        if not os.path.exists(os.path.split("./figs_pdf/" + name + ".pdf")[0]):
            os.makedirs(os.path.split("./figs_pdf/" + name + ".pdf")[0])
        fig.savefig("./figs_pdf/" + name + ".pdf")
    if png:
        if not os.path.exists(os.path.split("./figs_png/" + name + ".png")[0]):
            os.makedirs(os.path.split("./figs_png/" + name + ".png")[0])
        fig.savefig("./figs_png/" + name + ".png", dpi=dpi)


# sets the limits and ticks, adds grid and sets aspect ratio equal
def fixticks(ax, xmin, xmax, ymin, ymax, dx, dy, xstart=None, gridalpha=1):
    """sets the limits and ticks, adds grid and sets aspect ratio equal

    args:
        - ax to be used
        - xmin, xmax, ymin, ymax: log10-values of ticks
        - dx, dy: steps of ticks
        - xstart: min. x value for labelled ticks
        - gridalpha: opacity of grid
    """
    if xstart is None:
        xstart = xmin
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(10**np.arange(int(np.log10(xstart)), int(np.log10(xmax))+1, dx, dtype=float))

    ax.set_xticks(10**np.arange(int(np.log10(xmin)), int(np.log10(xmax))+1, 1, dtype=float), minor=True)
    ax.set_xticklabels(["" for i in np.arange(int(np.log10(xmin)), int(np.log10(xmax))+1, 1, dtype=float)], minor=True)
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(10**np.arange(int(np.log10(ymin)), int(np.log10(ymax))+1, dy, dtype=float))

    ax.set_yticks(10**np.arange(int(np.log10(ymin)), int(np.log10(ymax))+1, 1, dtype=float), minor=True)
    ax.set_yticklabels(["" for i in np.arange(int(np.log10(ymin)), int(np.log10(ymax))+1, 1, dtype=float)], minor=True)
    ax.grid(which="both", visible=True, alpha=gridalpha)
    ax.set_aspect("equal")
    return ax