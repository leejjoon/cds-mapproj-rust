import numpy as np
from mapproj.zenithal.sin import Sin

def test_sin_projection():
    sin = Sin()

    # Test projection
    xyz = (1.0, 0.0, 0.0)
    xy = sin.proj(xyz)
    assert np.allclose(xy, (0.0, 0.0))

    xyz = (0.5, 0.5, np.sqrt(0.5))
    xy = sin.proj(xyz)
    assert np.allclose(xy, (0.5, np.sqrt(0.5)))

    # Test projection of a point on the back hemisphere
    xyz = (-0.5, 0.5, 0.5)
    xy = sin.proj(xyz)
    assert np.isnan(xy).all()

    # Test unprojection
    xy = (0.0, 0.0)
    xyz = sin.unproj(xy)
    assert np.allclose(xyz, (1.0, 0.0, 0.0))

    xy = (0.5, 0.5)
    xyz = sin.unproj(xy)
    assert np.allclose(xyz, (np.sqrt(0.5), 0.5, 0.5))

    # Test unprojection of a point outside the projection plane
    xy = (1.1, 0.0)
    xyz = sin.unproj(xy)
    assert np.isnan(xyz).all()

    # Test projection and unprojection consistency
    xyz_original = (np.sqrt(0.5), 0.5, 0.5)
    xy = sin.proj(xyz_original)
    xyz_new = sin.unproj(xy)
    assert np.allclose(xyz_original, xyz_new)
