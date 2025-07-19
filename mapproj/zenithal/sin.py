import numpy as np
from typing import Tuple

class Sin:
    """Orthographic projection."""

    NAME = "Orthographic"
    WCS_NAME = "SIN"

    def __init__(self):
        pass

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        return -1.0, 1.0, -1.0, 1.0

    def proj(self, xyz: Tuple[float, float, float]) -> Tuple[float, float]:
        """Projects 3D coordinates to 2D."""
        if xyz[0] >= 0.0:
            return xyz[1], xyz[2]
        else:
            return np.nan, np.nan

    def unproj(self, xy: Tuple[float, float]) -> Tuple[float, float, float]:
        """Unprojects 2D coordinates to 3D."""
        r2 = xy[0]**2 + xy[1]**2
        if r2 <= 1.0:
            x = np.sqrt(1.0 - r2)
            return x, xy[0], xy[1]
        else:
            return np.nan, np.nan, np.nan

class SinSlant:
    """Slant Orthographic projection."""

    NAME = "Slant orthographic"
    WCS_NAME = "SIN"

    def __init__(self, xi: float, eta: float):
        self.xi = xi
        self.eta = eta
        self.tg2 = xi**2 + eta**2
        tmp = np.sqrt(1.0 + self.tg2)
        self.xp = -1.0 / tmp
        self.yp = -xi / tmp
        self.zp = -eta / tmp
        self.proj_bounds = (-1.0, 1.0 + xi * 2.0, -1.0, 1.0 + eta * 2.0)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        return self.proj_bounds

    def proj(self, xyz: Tuple[float, float, float]) -> Tuple[float, float]:
        """Projects 3D coordinates to 2D."""
        s = xyz[0] * self.xp + xyz[1] * self.yp + xyz[2] * self.zp
        if s <= 0.0:
            s = 1.0 - xyz[0]
            return xyz[1] + self.xi * s, xyz[2] + self.eta * s
        else:
            return np.nan, np.nan

    def unproj(self, pos: Tuple[float, float]) -> Tuple[float, float, float]:
        """Unprojects 2D coordinates to 3D."""
        x2d, y2d = pos
        r2 = x2d**2 + y2d**2
        s = self.xp + x2d * self.yp + y2d * self.yp
        v = (1.0 - self.xp * s)**2 + (x2d - self.yp * s)**2 + (y2d - self.zp * s)**2
        if v < 1.0:
            rp = self.xi * x2d + self.eta * y2d
            a = 1.0 + self.tg2
            b = 2.0 * (rp - self.tg2)
            c = r2 - 2.0 * rp + self.tg2 - 1.0
            x = (-b + np.sqrt(b**2 - 4.0 * a * c)) / (2.0 * a)
            tmp = 1.0 - x
            return x, x2d - self.xi * tmp, y2d - self.eta * tmp
        else:
            return np.nan, np.nan, np.nan
