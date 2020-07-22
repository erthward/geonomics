import tskit
import matplotlib.pyplot as pyplot

def ld_matrix_example(ts):
    ld_calc = tskit.LdCalculator(ts)
    A = ld_calc.r2_matrix()
    # Now plot this matrix.
    x = A.shape[0] / pyplot.rcParams["figure.dpi"]
    x = max(x, pyplot.rcParams["figure.figsize"][0])
    fig, ax = pyplot.subplots(figsize=(x, x))
    fig.tight_layout(pad=0)
    im = ax.imshow(A, interpolation="none", vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in "top", "bottom", "left", "right":
        ax.spines[s].set_visible(False)
    pyplot.gcf().colorbar(im, shrink=0.5, pad=0)
    #pyplot.savefig("ld.svg")
