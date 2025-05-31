import math
from mapproj import PyLonLat, PyImgXY, PySin, PyCenteredProjectionSin, PyWcsImgXY2ProjXY, PyImg2CelestialSinWcs # Changed import

def test_readme_example():
    # Define constants
    crpix1 = 382.00001513958
    crpix2 = 389.500015437603

    crval1_deg = 183.914583333
    crval2_deg = 36.3275

    cd11 = -2.7777777349544e-4
    cd22 = 2.77777773495436e-4

    # Set the projection
    sin_proj = PySin()
    # The Rust example creates CenteredProjection with Sin::default()
    # then sets the center.
    # PyCenteredProjectionSin takes PySin in constructor.
    # The original RsCenteredProjection takes the actual projection struct (e.g. RsSin), not a PySin.
    # The PyCenteredProjectionSin constructor takes PySin, and internally uses sin_proj.proj which is RsSin.
    centered_proj = PyCenteredProjectionSin(sin_proj)

    proj_center_lon_rad = math.radians(crval1_deg)
    proj_center_lat_rad = math.radians(crval2_deg)
    proj_center = PyLonLat(proj_center_lon_rad, proj_center_lat_rad)

    centered_proj.set_proj_center_from_lonlat(proj_center)

    # In Rust: ImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, 0.0, 0.0, cd22);
    # PyWcsImgXY2ProjXY constructor is from_cd
    wcs_transform = PyWcsImgXY2ProjXY(crpix1, crpix2, cd11, 0.0, 0.0, cd22)

    # In Rust: Img2Celestial::new(img2proj, proj);
    # PyImg2CelestialSinWcs constructor takes PyCenteredProjectionSin and PyWcsImgXY2ProjXY
    img2celestial = PyImg2CelestialSinWcs(centered_proj, wcs_transform)
    # We could have set the projection center here instead of previously:
    # img2celestial.set_proj_center_from_lonlat(proj_center)


    # Use to project, unproject coordinates:
    # - we choose on purpose position in the image of the projection center
    img_coo_input_py = PyImgXY(crpix1, crpix2)

    lonlat_py_opt = img2celestial.img2lonlat(img_coo_input_py)
    assert lonlat_py_opt is not None
    lonlat_py = lonlat_py_opt

    assert math.isclose(lonlat_py.lon, proj_center.lon, abs_tol=1e-14)
    assert math.isclose(lonlat_py.lat, proj_center.lat, abs_tol=1e-14)

    # The variable name in Rust example is img_coo_input, but it's a result of lonlat2img
    # Renaming to avoid confusion
    img_coo_output_py_opt = img2celestial.lonlat2img(lonlat_py)
    assert img_coo_output_py_opt is not None
    img_coo_output_py = img_coo_output_py_opt

    # The original Rust test compares img_coo_input.x() with img_coo_input.x()
    # This is likely a typo and should be img_coo_output.x() with img_coo_input.x()
    assert math.isclose(img_coo_output_py.x, img_coo_input_py.x, abs_tol=1e-14)
    assert math.isclose(img_coo_output_py.y, img_coo_input_py.y, abs_tol=1e-14)

def test_sin_out_of_bounds():
    # Define constants for a simple case
    crpix1 = 0.0
    crpix2 = 0.0
    # Projection center at (0, 0)
    crval1_deg = 0.0
    crval2_deg = 0.0
    # Simple scaling: 1 deg per pixel, no rotation
    # cd11 = 1.0 means 1 deg/pix on X, cd22 = 1.0 means 1 deg/pix on Y
    # However, the projection plane for SIN is in range [-1, 1]
    # We need to define CD matrix that maps pixel coordinates to something meaningful in [-1,1]
    # Let's assume a 100x100 pixel image covers the entire front hemisphere.
    # So pixel (50,50) is (0,0) in projection plane, pixel (0,50) is (-1,0)
    # cd11 will map pixel span to projection plane span.
    # If image is 100 pixels wide, and proj plane is 2 units wide (-1 to 1), then 1 pixel = 0.02 proj units.
    cd11 = 0.02
    cd22 = 0.02

    sin_proj = PySin()
    centered_proj = PyCenteredProjectionSin(sin_proj)

    proj_center_lon_rad = math.radians(crval1_deg)
    proj_center_lat_rad = math.radians(crval2_deg)
    proj_center = PyLonLat(proj_center_lon_rad, proj_center_lat_rad)

    centered_proj.set_proj_center_from_lonlat(proj_center)

    wcs_transform = PyWcsImgXY2ProjXY(crpix1, crpix2, cd11, 0.0, 0.0, cd22)
    img2celestial = PyImg2CelestialSinWcs(centered_proj, wcs_transform)

    # Point on the back of the sphere (opposite to projection center 0,0)
    # e.g., Lon = 180 deg (pi rad), Lat = 0 deg
    lon_back_rad = math.pi
    lat_back_rad = 0.0
    lonlat_back = PyLonLat(lon_back_rad, lat_back_rad)

    img_coo_opt = img2celestial.lonlat2img(lonlat_back)
    assert img_coo_opt is None

    # Test img2lonlat with coordinates that would be far out if projected
    # For SIN, the projection plane is bounded by x,y in [-1, 1].
    # If cd11 = 0.02, then an image x coordinate of 100 would map to:
    # (100 - crpix1) * cd11 = (100 - 0) * 0.02 = 2.0 in projection plane x.
    # This is outside the SIN valid area.
    img_far_out = PyImgXY(100.0, 0.0) # x_proj = 2.0
    lonlat_opt = img2celestial.img2lonlat(img_far_out)
    assert lonlat_opt is None

    img_far_out_y = PyImgXY(0.0, 100.0) # y_proj = 2.0
    lonlat_opt_y = img2celestial.img2lonlat(img_far_out_y)
    assert lonlat_opt_y is None

    # A point that should be valid
    # (0,0) in image coords, with crpix1=0, crpix2=0, should be (0,0) in proj plane
    # which is the projection center, (0,0) lon/lat
    img_center = PyImgXY(0.0, 0.0)
    lonlat_center_opt = img2celestial.img2lonlat(img_center)
    assert lonlat_center_opt is not None
    assert math.isclose(lonlat_center_opt.lon, 0.0, abs_tol=1e-9)
    assert math.isclose(lonlat_center_opt.lat, 0.0, abs_tol=1e-9)
