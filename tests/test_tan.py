import numpy as np
from mapproj.zenithal.tan import Tan

def test_tan_projection():
    tan = Tan()

    # Test projection
    xyz = (1.0, 0.0, 0.0)
    xy = tan.proj(xyz)
    assert np.allclose(xy, (0.0, 0.0))

    xyz = (0.5, 0.5, 0.5)
    xy = tan.proj(xyz)
    assert np.allclose(xy, (1.0, 1.0))

    # Test projection of a point on the back hemisphere
    xyz = (-0.5, 0.5, 0.5)
    xy = tan.proj(xyz)
    assert np.isnan(xy).all()

    # Test unprojection
    xy = (0.0, 0.0)
    xyz = tan.unproj(xy)
    assert np.allclose(xyz, (1.0, 0.0, 0.0))

    xy = (1.0, 1.0)
    xyz = tan.unproj(xy)
    expected = 1.0 / np.sqrt(3)
    assert np.allclose(xyz, (expected, expected, expected))

    # Test projection and unprojection consistency
    xyz_original = (1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3))
    xy = tan.proj(xyz_original)
    xyz_new = tan.unproj(xy)
    assert np.allclose(xyz_original, xyz_new)
