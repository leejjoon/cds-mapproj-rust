from astropy.io import fits
from astropy.wcs import WCS

# example header from https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
s = """CTYPE1  = 'RA---TAN-SIP' / RA---TAN with distortion
CTYPE2  = 'DEC--TAN-SIP' / DEC--TAN with distortion
CRVAL1  = 248.588845714283 / [deg] RA at CRPIX1,CRPIX2
CRVAL2  = -46.7664545714299 / [deg] DEC at CRPIX1,CRPIX2
CD1_1   = 3.79787003748175E-05
CD1_2   = 0.000336011759341397
CD2_1   = 0.000337467809645458
CD2_2   = -3.6366440945003E-05
NAXIS   = 2
NAXIS1  = 256
NAXIS2  = 256
CRPIX1  = 128.
CRPIX2  = 128.
A_ORDER = 2 / polynomial order, axis 1
A_0_2   = 1.569E-05 / distortion coefficient
A_1_1   = 5.232E-05 / distortion coefficient
A_2_0   = 3.31E-05 / distortion coefficient
B_ORDER = 2 / polynomial order, axis 2
B_0_2   = 4.172E-05 / distortion coefficient
B_1_1   = 2.213E-05 / distortion coefficient
B_2_0   = -9.819E-07 / distortion coefficient
AP_ORDER= 2 / polynomial order, axis 1
AP_0_1  = 5.677E-05 / distortion coefficient
AP_0_2  = -1.569E-05 / distortion coefficient
AP_1_0  = 5.871E-05 / distortion coefficient
AP_1_1  = -5.231E-05 / distortion coefficient
AP_2_0  = -3.309E-05 / distortion coefficient
BP_ORDER= 2 / polynomial order, axis 2
BP_0_1  = 4.432E-05 / distortion coefficient
BP_0_2  = -4.172E-05 / distortion coefficient
BP_1_0  = 2.091E-05 / distortion coefficient
BP_1_1  = -2.213E-05 / distortion coefficient
BP_2_0  = 9.814E-07 / distortion coefficient
"""
import re
p_sip = re.compile(r"[AB](P?)\_(\w)+")
# p_sip = re.compile(r"[AB](P?)\_\d_\d")
# p_sip_order = re.compile(r"[AB](P?)\_ORDER")

ss = [s1 for s1 in s.split("\n")]
# cards = [fits.Card.fromstring(s1).verify("fix") for s1 in s.split("\n")]
cards = [fits.Card.fromstring(s1) for s1 in s.split("\n")]

h = fits.Header(cards)
w = WCS(h)

if False:
    cards_no_sip = [c for c in cards if not p_sip.match(c.keyword)]
    h = fits.Header(cards_no_sip)
    h["CTYPE1"] = h["CTYPE1"].replace("-SIP", "")
    h["CTYPE2"] = h["CTYPE2"].replace("-SIP", "")
    w_no_sip = WCS(h)

import numpy as np
# m = np.mgrid[-64:64+1:64, -32:32+1:32]

def test1():
    x0, y0 = np.meshgrid(np.linspace(-64, 64, 3), np.linspace(-32, 32, 3))
    x = h["CRPIX1"] + x0
    y = h["CRPIX2"] + y0

    u, v = w.sip_pix2foc(x, y, 1)
    xx,yy = w.sip_foc2pix(u, v, 1)

    xx0 = x0 - h["CRPIX1"]
    yy0 = y0 - h["CRPIX2"]

def test2():
    # from pix to world

    crval1 = h["CRVAL1"]
    crval2 = h["CRVAL2"]
    crpix1 = h["CRPIX1"]
    crpix2 = h["CRPIX2"]

    cd = 4.e-5

    u_range = np.linspace(-64, 64, 3)
    v_range = np.linspace(-32, 32, 3)
    u, v = np.meshgrid(u_range, v_range)

    x = crpix1 + u
    y = crpix2 + v

    uf, vg = w.sip_pix2foc(x, y, 1) # uf = u + f, vg = v + g

    print("// from pixel to foc")
    print(f"let u_range = {u_range[0]}..={u_range[-1]}")
    print(f"let v_range = {v_range[0]}..={v_range[-1]}")
    for l in get_rust_vec(h):
        print(l)

    for i, (u1, v1, f1, g1) in enumerate(zip(u.flat, v.flat, uf.flat, vg.flat)):
        t = ",".join(map(str, [u1, v1, f1, g1]))
        print(f"let uvfg = ({t});")

    lon0, lat0 = w.wcs_pix2world(uf + crpix1, vg + crpix2, 1)

    lon1, lat1 = w.all_pix2world(x, y, 1)

    assert np.allclose([lon0, lat0], [lon1, lat1], atol=1.e-3*cd) # in radian


def test3():
    # from world to pix

    crval1 = h["CRVAL1"]
    crval2 = h["CRVAL2"]
    crpix1 = h["CRPIX1"]
    crpix2 = h["CRPIX2"]

    cd = 4.e-4
    dlon, dlat = np.meshgrid(cd * np.linspace(-64, 64, 3),
                             cd * np.linspace(-32, 32, 3))
    lon = crval1 + dlon
    lat = crval2 + dlat

    uF, vG = w.wcs_world2pix(lon, lat, 1)

    x0, y0 = w.sip_foc2pix(uF - crpix1, vG - crpix2, 1)

    x1, y1 = w.all_world2pix(lon, lat, 1)

    assert np.allclose([x0, y0], [x1, y1], atol=1.e-2) # in pixel

def get_rust_vec(h):
    # "let a_coeff = SipCoeff::new(Box::new([0., 0., 1.]));"

    rust_codes = []
    for ab in ["A", "B", "AP", "BP"]:
        order = h[f"{ab}_ORDER"]
        c = []
        for p in range(order + 1):
            for q in range(order - p + 1):
                k = f"{ab}_{p}_{q}"
                c.append(h.get(k, 0.))
                # print(p, q, h.get(k, 0))

        ab_coeff = ", ".join(map(str, c))
        rust_code = f"let {ab.lower()}_coeff = SipCoeff::new(Box::new([{ab_coeff}]));"

        rust_codes.append(rust_code)

    return rust_codes
