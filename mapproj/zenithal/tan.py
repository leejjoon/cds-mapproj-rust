import numpy as np
from typing import Tuple

class Tan:
    """Gnomonic projection."""

    NAME = "Gnomonic"
    WCS_NAME = "TAN"

    def __init__(self):
        pass

    @property
    def bounds(self) -> None:
        return None

    def proj(self, xyz: Tuple[float, float, float]) -> Tuple[float, float]:
        """Projects 3D coordinates to 2D."""
        if xyz[0] > 0.0:
            return xyz[1] / xyz[0], xyz[2] / xyz[0]
        else:
            return np.nan, np.nan

    def unproj(self, pos: Tuple[float, float]) -> Tuple[float, float, float]:
        """Unprojects 2D coordinates to 3D."""
        x = 1.0 / np.sqrt(1.0 + pos[0]**2 + pos[1]**2)
        return x, pos[0] * x, pos[1] * x
