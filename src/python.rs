use pyo3::prelude::*;
// wrap_pyfunction seems unused based on previous warnings, let's remove it for now.
// use pyo3::wrap_pyfunction;

// Items from lib.rs are now directly available via crate::
use crate::{LonLat as RsLonLat, ImgXY as RsImgXY};

#[pyclass(name = "LonLat")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyLonLat {
  #[pyo3(get, set)]
  lon: f64,
  #[pyo3(get, set)]
  lat: f64,
}

#[pymethods]
impl PyLonLat {
  #[new]
  pub fn new(lon: f64, lat: f64) -> Self {
    Self {lon, lat}
  }
}

impl From<RsLonLat> for PyLonLat {
  fn from(lonlat: RsLonLat) -> Self {
    Self {
      lon: lonlat.lon(),
      lat: lonlat.lat(),
    }
  }
}

impl From<PyLonLat> for RsLonLat {
  fn from(lonlat: PyLonLat) -> Self {
    RsLonLat::new(lonlat.lon, lonlat.lat)
  }
}

#[pyclass(name = "ImgXY")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyImgXY {
  #[pyo3(get, set)]
  x: f64,
  #[pyo3(get, set)]
  y: f64,
}

#[pymethods]
impl PyImgXY {
  #[new]
  pub fn new(x: f64, y: f64) -> Self {
    Self { x, y }
  }
}

impl From<RsImgXY> for PyImgXY {
  fn from(img_xy: RsImgXY) -> Self {
    Self {
      x: img_xy.x(),
      y: img_xy.y(),
    }
  }
}

impl From<PyImgXY> for RsImgXY {
  fn from(img_xy: PyImgXY) -> Self {
    RsImgXY::new(img_xy.x, img_xy.y)
  }
}

use crate::zenithal::sin::Sin as RsSin; // Corrected path

#[pyclass(name = "Sin")]
#[derive(Clone)]
pub struct PySin {
  proj: RsSin,
}

#[pymethods]
impl PySin {
  #[new]
  pub fn default() -> Self {
    Self { proj: RsSin::default() }
  }
}

// CenteredProjection is used. Projection was unused.
use crate::CenteredProjection as RsCenteredProjection;
// use crate::Projection; // Unused import


#[pyclass(name = "CenteredProjectionSin")]
#[derive(Clone)]
pub struct PyCenteredProjectionSin {
  proj: RsCenteredProjection<RsSin>,
}

#[pymethods]
impl PyCenteredProjectionSin {
  #[new]
  pub fn new(sin_proj: PySin) -> Self {
    Self { proj: RsCenteredProjection::new(sin_proj.proj) }
  }

  pub fn set_proj_center_from_lonlat(&mut self, lonlat: PyLonLat) {
    self.proj.set_proj_center_from_lonlat(&RsLonLat::from(lonlat));
  }
}

use crate::img2proj::WcsImgXY2ProjXY as RsWcsImgXY2ProjXY;

#[pyclass(name = "WcsImgXY2ProjXY")]
#[derive(Clone)]
pub struct PyWcsImgXY2ProjXY {
  transform: RsWcsImgXY2ProjXY,
}

#[pymethods]
impl PyWcsImgXY2ProjXY {
  #[new]
  pub fn from_cd(
    crpix1: f64, crpix2: f64,
    cd11: f64, cd12: f64,
    cd21: f64, cd22: f64,
  ) -> Self {
    Self { transform: RsWcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22) }
  }
}

use crate::img2celestial::Img2Celestial as RsImg2Celestial;
// use crate::img2proj::ImgXY2ProjXY as RsImgXY2ProjXY; // Unused import

#[pyclass(name = "Img2CelestialSinWcs")]
pub struct PyImg2CelestialSinWcs {
  img2celestial: RsImg2Celestial<RsSin, RsWcsImgXY2ProjXY>,
}

#[pymethods]
impl PyImg2CelestialSinWcs {
  #[new]
  pub fn new(centered_proj: PyCenteredProjectionSin, wcs_transform: PyWcsImgXY2ProjXY) -> Self {
    Self {
      img2celestial: RsImg2Celestial::new(wcs_transform.transform, centered_proj.proj)
    }
  }

  pub fn set_proj_center_from_lonlat(&mut self, lonlat: PyLonLat) {
    self.img2celestial.set_proj_center_from_lonlat(&RsLonLat::from(lonlat));
  }

  pub fn lonlat2img(&self, lonlat: PyLonLat) -> PyResult<Option<PyImgXY>> {
    match self.img2celestial.lonlat2img(&RsLonLat::from(lonlat)) {
      Some(img_xy) => Ok(Some(PyImgXY::from(img_xy))),
      None => Ok(None),
    }
  }

  pub fn img2lonlat(&self, img_xy: PyImgXY) -> PyResult<Option<PyLonLat>> {
    match self.img2celestial.img2lonlat(&RsImgXY::from(img_xy)) {
      Some(lonlat) => Ok(Some(PyLonLat::from(lonlat))),
      None => Ok(None),
    }
  }
}

#[pymodule]
fn mapproj(_py: Python, m: &PyModule) -> PyResult<()> { // Changed function name to mapproj
    m.add_class::<PyLonLat>()?;
    m.add_class::<PyImgXY>()?;
    m.add_class::<PySin>()?;
    m.add_class::<PyCenteredProjectionSin>()?;
    m.add_class::<PyWcsImgXY2ProjXY>()?;
    m.add_class::<PyImg2CelestialSinWcs>()?;
    Ok(())
}
