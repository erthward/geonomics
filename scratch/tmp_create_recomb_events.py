import numpy as np
import matplotlib.pyplot as plt


def _create_recomb_events(n_events, positions, alpha=None, beta=None,
                          r=None):
    """
    NOTE: If positions and/or r are provided,
    they must already be sorted!
    """
    positions = [*positions]
    positions = np.array(positions)
    events = []
    if r is not None:
        assert len(r) == len(positions), ("Lengths of provided recombination "
                                          "rates and recombination "
                                          "positions don't match!")
    elif alpha is not None and beta is not None:
        r = np.clip(np.random.beta(a=alpha, b=beta, size=len(positions)),
                    a_min=0, a_max=0.5)
    elif alpha is not None:
        r = np.ones(len(positions)) * alpha
    else:
        r = np.ones(len(positions)) * 0.5
    events = [positions[np.where(np.random.binomial(1, r))] for _ in range(
                                                                    n_events)]
    return (np.array([*positions]), r, events)
