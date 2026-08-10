"""Numerically stable statistical identities shared by harmonisation steps."""

import math


def two_sided_negative_log10_p_from_z(values):
    """Return ``-log10(2 * P(N(0,1) > |Z|))`` without tail underflow."""
    import numpy as np
    from scipy.special import log_ndtr

    z_values = np.asarray(values, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return -(
            math.log(2.0) + log_ndtr(-np.abs(z_values))
        ) / math.log(10.0)


__all__ = ["two_sided_negative_log10_p_from_z"]
