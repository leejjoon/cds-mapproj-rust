//! Module containing the structure to convert back on forth 
//! from Image coordinates to Intermediate coordinates, i.e. coordinates in the projection plane.

use std::ops::RangeInclusive;
use crate::{
  ImgXY, ProjXY,
  sip::Sip
};

/// Transform the XY coordinates in the projection plane in a pixel coordinates in an image.
pub trait ProjXY2ImgXY {
  
  /// Transforms intermediate world coordinates -- i.e. coordinates in the canonical projection 
  /// plane -- to pixel coordinates in an image.
  /// # Params
  /// * `xy` coordinates in the canonical projection plane.
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY>;
}

/// Transform the pixel coordinates in an image to the XY coordinates in the projection plane.
pub trait ImgXY2ProjXY {
  
  /// Type of the inverse transformation.
  type T: ProjXY2ImgXY;

  /// Transforms pixel coordinates to the intermediate world coordinates, i.e.
  ///  to coordinates in the canonical projection plane.
  /// # Params
  /// * `xy`: pixel coordinates (no units)
  fn img2proj(&self, xy: &ImgXY) -> ProjXY;
  
  /// Provide the inverse transformation.
  fn inverse(&self) -> Self::T;
}

#[derive(Clone)]
pub struct BasicImgXY2ProjXY {
  /// Projection x-axis coordinate at the image X-axis center.
  center_px: f64,
  /// Projection y-axis coordinate at the image Y-axis center.
  center_py: f64,
  /// Center of the image in the X-axis (= half the number of pixels along the X-axis).
  center_x: f64,
  /// Center of the image in the Y-axis (= half the number of pixels along the Y-axis).
  center_y: f64,
  /// Size of a pixel along the projection plane the x-axis.
  scale_x: f64,
  /// Size of a pixel along the projection plane the y-axis.
  scale_y: f64,
}

impl BasicImgXY2ProjXY {
  
  /// We assume the origin to be on the bottom left corner of the image,
  /// with the `x` increasing from left to right and `y` increasing bottom-up.
  /// # Remark
  /// * for a PNG, the origin is the top left corner and `y` is increasing top-down.
  /// * for the east to be towards the left of an image, the `x` axis as to be reversed.
  /// # Params
  /// * `img_size`: `(size_x, size_y)` number of pixels in each dimension
  /// * `proj_bounds`: `(bounds_x, bounds_y)` boundaries of the projection domain
  pub fn from(
    img_size: (u16, u16),
    proj_bounds: (&RangeInclusive<f64>, &RangeInclusive<f64>)
  ) -> Self {
    // x-axis
    let img_size_x = img_size.0 as f64;
    let center_x = 0.5 * (img_size_x - 1.0);
    let xp_min = proj_bounds.0.start();
    let xp_max = proj_bounds.0.end();
    let scale_x = (xp_max - xp_min) / img_size_x as f64;
    let center_px = 0.5 * (xp_max + xp_min);
    // y-axis
    let img_size_y = img_size.1 as f64;
    let center_y = 0.5 * (img_size_y - 1.0);
    let yp_min = proj_bounds.1.start();
    let yp_max = proj_bounds.1.end();
    let scale_y = (yp_max - yp_min) / img_size_y as f64;
    let center_py = 0.5 * (yp_max + yp_min);
    // Create
    Self {
      center_px,
      center_py,
      center_x,
      center_y,
      scale_x,
      scale_y
    }
  }
}

impl ImgXY2ProjXY for BasicImgXY2ProjXY {
  type T = Self;

  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    let proj_x = self.center_px + (xy.x as f64 - self.center_x) * self.scale_x;
    let proj_y = self.center_py + (xy.y as f64 - self.center_y) * self.scale_y;
    ProjXY::new(proj_x, proj_y)
  }

  fn inverse(&self) ->  Self::T {
    self.clone()
  }
}

impl ProjXY2ImgXY for BasicImgXY2ProjXY {

  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.center_x + (xy.x - self.center_px) / self.scale_x;
    let y = self.center_y + (xy.y - self.center_py) / self.scale_y;
    Some(ImgXY::new(x, y))
  }
}


/// Specific implementation for PNG (top left origin) with East towards the left.
#[derive(Clone)]
pub struct ReversedEastPngImgXY2ProjXY{
  y_img_max: f64,
  b: BasicImgXY2ProjXY,
}

impl ReversedEastPngImgXY2ProjXY {
  pub fn from(
    img_size: (u16, u16),
    proj_bounds: (&RangeInclusive<f64>, &RangeInclusive<f64>)
  ) -> Self {
    Self{
      y_img_max: (img_size.1 - 1) as f64,
      b: BasicImgXY2ProjXY::from(img_size, proj_bounds) 
    }
  }
}


impl ImgXY2ProjXY for ReversedEastPngImgXY2ProjXY {
  type T = Self;

  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    let proj_x = self.b.center_px + (xy.x - self.b.center_x) * self.b.scale_x;
    let proj_y = self.b.center_py + ((self.y_img_max - xy.y) - self.b.center_y) * self.b.scale_y;
    ProjXY::new(-proj_x, proj_y)
  }

  fn inverse(&self) ->  Self::T {
    self.clone()
  }
}

impl ProjXY2ImgXY for ReversedEastPngImgXY2ProjXY {

  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.b.center_x + (-xy.x - self.b.center_px) / self.b.scale_x;
    let y = self.b.center_y + (xy.y - self.b.center_py) / self.b.scale_y;
    Some(ImgXY::new(x, self.y_img_max - y))
  }
}



/// Struct allowing to transform the pixel coordinates in an image to the XY coordinates 
/// in the projection plane.
/// The three constructors are each associated with one of the three convention
/// describe in the FITS paper: CDij, CDELTi + PCij, CDELTi + CROTA2.
#[derive(Clone)]
pub struct WcsImgXY2ProjXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64, crpix2: f64,
  /// Rotation (no units) combined with a scale (in radians) matrix.
  cd11: f64, cd12: f64,
  cd21: f64, cd22: f64,
}

impl WcsImgXY2ProjXY {
  
  /// Create a struct from the `CDij` convention.
  /// # Params
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `cd11`: value of the `CD11` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  /// * `cd12`: value of the `CD12` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd21`: value of the `CD21` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd22`: value of the `CD22` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  pub fn from_cd(
    crpix1: f64, crpix2: f64,
    cd11: f64, cd12: f64, 
    cd21: f64, cd22: f64,
  ) -> Self {
    eprintln!("[mapproj DEBUG WcsImgXY2ProjXY::from_cd] crpix1: {:.15}, crpix2: {:.15}", crpix1, crpix2);
    Self {
      crpix1, crpix2,
      cd11: cd11.to_radians(),
      cd12: cd12.to_radians(),
      cd21: cd21.to_radians(),
      cd22: cd22.to_radians(),
    }
  }

  /// Create a struct from the `CDELTi + PCij` convention.
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `pc11`: value of the `PC11` keyword (element of a rotation matrix, no units), use 1 as default value
  /// * `pc12`: value of the `PC12` keyword (element of a rotation matrix, no units), use 0 as default value
  /// * `pc21`: value of the `PC21` keyword (element of a rotation matrix, no units), use 0 as default value
  /// * `pc22`: value of the `PC22` keyword (element of a rotation matrix, no units), use 1 as default value
  /// * `cdelt1`: value of the `CDELT1` keyword, in degrees
  /// * `cdelt2`: value of the `CDELT2` keyword, in degrees
  pub fn from_pc(
    crpix1: f64, crpix2: f64,
    pc11: f64,  pc12: f64,
    pc21: f64,  pc22: f64,
    cdelt1: f64, cdelt2: f64
  ) -> Self {
    eprintln!("[mapproj DEBUG WcsImgXY2ProjXY::from_pc] crpix1: {:.15}, crpix2: {:.15}", crpix1, crpix2);
    Self::from_cd(
      crpix1, crpix2,
      cdelt1 * pc11, cdelt1 * pc12,
      cdelt2 * pc21, cdelt2 * pc22,
    )
  }

  /// Create a struct from the `CDELTi + CROTA2` convention.
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `crota2`: value of the `CROTA2` keyword, in degrees
  /// * `cdelt1`: value of the `CDELT1` keyword, in degrees
  /// * `cdelt2`: value of the `CDELT2` keyword, in degrees
  pub fn from_cr(
    crpix1: f64, crpix2: f64,
    crota2: f64,
    cdelt1: f64, cdelt2: f64
  ) -> Self {
    eprintln!("[mapproj DEBUG WcsImgXY2ProjXY::from_cr] crpix1: {:.15}, crpix2: {:.15}", crpix1, crpix2);
    let (sinc, cosc) = crota2.to_radians().sin_cos();
    Self::from_cd(
      crpix1, crpix2, 
      cdelt1 * cosc, cdelt1 * sinc,
      -cdelt2 * sinc, cdelt2 * cosc
    )
  }
  
}

impl ImgXY2ProjXY for WcsImgXY2ProjXY {
  
  type T = WcsProjXY2ImgXY;
  
  /// Transform the pixel coordinates to the intermediate world coordinates
  /// (or native spherical coordinates) by applying first a translation
  /// (given the `CRPIXi` keywords value) and then a rotation plus a scale
  /// (given the `CDij` keywords values).
  /// # Params
  /// * `imgXY`: pixel coordinates (no units) to be transformed intermediate world coordinates
  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    eprintln!("[mapproj DEBUG WcsImgXY2ProjXY::img2proj] self.crpix1: {:.15}, self.crpix2: {:.15}", self.crpix1, self.crpix2);
    eprintln!("[mapproj DEBUG WcsImgXY2ProjXY::img2proj] input xy.x: {:.15}, xy.y: {:.15}", xy.x, xy.y);
    // Translation
    let x = xy.x - self.crpix1;
    let y = xy.y - self.crpix2;
    // Rotation + scale
    ProjXY::new(
      self.cd11 * x + self.cd12 * y,
      self.cd21 * x + self.cd22 * y
    )
  }
  
  
  fn inverse(&self) -> Self::T {
    // Compute the determinant of the CD matrix
    let det = self.cd11 * self.cd22 - self.cd12 * self.cd21;
    // Compute the coefficient of the inverse matrix
    WcsProjXY2ImgXY {
      crpix1: self.crpix1,
      crpix2: self.crpix2,
      icd11:  self.cd22 / det,
      icd12: -self.cd21 / det,
      icd21: -self.cd12 / det,
      icd22:  self.cd11 / det,
    }
  }
}

/// Struct allowing to transform the pixel coordinates in an image to the XY coordinates 
/// in the projection plane.
/// The three constructors are each associated with one of the three convention
/// describe in the FITS paper: CDij, CDELTi + PCij, CDELTi + CROTA2.
pub struct WcsWithSipImgXY2ProjXY {
  /// Regular transformation
  wcs: WcsImgXY2ProjXY,
  /// SIP transformation
  sip: Sip,
}

impl WcsWithSipImgXY2ProjXY {
  
  /// Add SIP convention to a regular WCS transformation.
  pub fn new(wcs: WcsImgXY2ProjXY, sip: Sip) -> Self {
    eprintln!("[mapproj DEBUG WcsWithSipImgXY2ProjXY::new] wcs.crpix1: {:.15}, wcs.crpix2: {:.15}", wcs.crpix1, wcs.crpix2);
    Self { wcs, sip }
  }
  
}

impl ImgXY2ProjXY for WcsWithSipImgXY2ProjXY {

  type T = WcsWithSipProjXY2ImgXY;

  /// Transform the pixel coordinates to the intermediate world coordinates
  /// (or native spherical coordinates) by applying first a translation
  /// (given the `CRPIXi` keywords value), then a rotation plus a scale
  /// (given the `CDij` keywords values), and finally the SIP distortion.
  /// # Params
  /// * `imgXY`: pixel coordinates (no units) to be transformed to intermediate world coordinates
  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    eprintln!("[mapproj DEBUG WcsWithSipImgXY2ProjXY::img2proj] self.wcs.crpix1: {:.15}, self.wcs.crpix2: {:.15}", self.wcs.crpix1, self.wcs.crpix2);
    eprintln!("[mapproj DEBUG WcsWithSipImgXY2ProjXY::img2proj] input xy.x: {:.15}, xy.y: {:.15}", xy.x, xy.y);
    // Compute 1-indexed pixel offsets from CRPIX (FITS header value).
    // SIP polynomials f(u,v), g(u,v), AP(u,v), BP(u,v) take u = x_image - CRPIX1 
    // and v = y_image - CRPIX2 as arguments, where x_image, y_image are 
    // 1-based pixel coordinates and CRPIX1, CRPIX2 are 1-based from FITS.
    // The input `xy` from wcs-rs is 0-based.
    // So, u_sip_arg = (xy.x + 1.0) - self.wcs.crpix1
    // And v_sip_arg = (xy.y + 1.0) - self.wcs.crpix2
    // These are the 1-indexed offsets required by SIP polynomials and CD matrix.
    let u_sip_arg = (xy.x + 0.0) - self.wcs.crpix1;
    let v_sip_arg = (xy.y + 0.0) - self.wcs.crpix2;
    eprintln!("[mapproj SIP INV DEBUG] u_sip_arg: {:.15e}, v_sip_arg: {:.15e}", u_sip_arg, v_sip_arg);

    // Apply inverse SIP distortion (AP and BP polynomials) using the new public methods.
    // The arguments to ap/bp are the same 1-indexed offsets u_sip_arg, v_sip_arg.
    let ap_val = self.sip.ap(u_sip_arg, v_sip_arg).unwrap_or(0.0);
    let bp_val = self.sip.bp(u_sip_arg, v_sip_arg).unwrap_or(0.0);
    eprintln!("[mapproj SIP INV DEBUG] ap_val: {:.15e}, bp_val: {:.15e}", ap_val, bp_val);

    // Calculate corrected pixel offsets u' and v' (FITS WCS Paper I, Eq. 6)
    let u_corrected = u_sip_arg + ap_val;
    let v_corrected = v_sip_arg + bp_val;
    eprintln!("[mapproj SIP INV DEBUG] u_corrected: {:.15e}, v_corrected: {:.15e}", u_corrected, v_corrected);

    // Linear transformation part: apply CD matrix to corrected pixel offsets
    let u_final_iwc = self.wcs.cd11 * u_corrected + self.wcs.cd12 * v_corrected;
    let v_final_iwc = self.wcs.cd21 * u_corrected + self.wcs.cd22 * v_corrected;
    eprintln!("[mapproj SIP INV DEBUG] u_final_iwc (CD*u'): {:.15e}, v_final_iwc (CD*v'): {:.15e}", u_final_iwc, v_final_iwc);

    ProjXY::new(u_final_iwc, v_final_iwc)
  }
  
  fn inverse(&self) -> Self::T {
    // Compute the inverse CD matrix
    let det = self.wcs.cd11 * self.wcs.cd22 - self.wcs.cd12 * self.wcs.cd21;
    // Compute the coefficient of the inverse matrix
    WcsWithSipProjXY2ImgXY {
      // The inverse WCS transformation (includes inverse CD matrix and CRPIX offset)
      wcs: WcsProjXY2ImgXY {
        crpix1: self.wcs.crpix1,
        crpix2: self.wcs.crpix2,
        // Inverse of CD matrix
        icd11: self.wcs.cd22 / det,
        icd12: -self.wcs.cd21 / det,
        icd21: -self.wcs.cd12 / det,
        icd22: self.wcs.cd11 / det,
      },
      // The same SIP coefficients are used for both forward and inverse transformations
      sip: self.sip.clone()
    }
  }
}



pub struct WcsProjXY2ImgXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64, crpix2: f64,
  /// Rotation (no units) combined with a scale (in radians) matrix.
  icd11: f64, icd12: f64,
  icd21: f64, icd22: f64,
}

impl ProjXY2ImgXY for WcsProjXY2ImgXY {
  
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    eprintln!("[mapproj DEBUG WcsProjXY2ImgXY::proj2img] self.crpix1: {:.15}, self.crpix2: {:.15}", self.crpix1, self.crpix2);
    eprintln!("[mapproj DEBUG WcsProjXY2ImgXY::proj2img] input xy.x: {:.15}, xy.y: {:.15}", xy.x, xy.y);
    let x = self.icd11 * xy.x + self.icd12 * xy.y + self.crpix1;
    let y = self.icd21 * xy.x + self.icd22 * xy.y + self.crpix2;
    Some(ImgXY::new(x, y))
   }
}


pub struct WcsWithSipProjXY2ImgXY {
  /// Regular transformation
  wcs: WcsProjXY2ImgXY,
  /// SIP transformation
  sip: Sip,
}

impl ProjXY2ImgXY for WcsWithSipProjXY2ImgXY {

  /// Transform from intermediate world coordinates to pixel coordinates by first
  /// removing the SIP distortion in world coordinates, then applying the inverse
  /// CD matrix, and finally adding the CRPIX offset.
  /// # Params
  /// * `xy`: intermediate world coordinates to be transformed to pixel coordinates
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    // 1. Perform the linear part of the transformation (inverse CD matrix + CRPIX offset)
    // This gives absolute pixel coordinates as if there were no SIP distortion.
    let img_xy_linear_absolute = self.wcs.proj2img(xy)?;

    // 2. Convert these absolute linear pixel coordinates to offsets from CRPIX.
    // The self.sip.inverse(u, v) method (which applies AP and BP polynomials)
    // Now, img_xy_linear contains pixel coordinates *as if there were no SIP distortion*.
    // These are effectively the (u, v) + CRPIX for the SIP polynomials A and B.
    // We need to calculate u and v (offsets from CRPIX)
    let u = img_xy_linear_absolute.x - self.wcs.crpix1; // This is u_intermediate_scaled
    let v = img_xy_linear_absolute.y - self.wcs.crpix2; // This is v_intermediate_scaled

    // Apply forward SIP distortion A and B
    // x_final = u + CRPIX + A(u,v)
    // y_final = v + CRPIX + B(u,v)
    // Since img_xy_linear_absolute.x() = u + CRPIX1_0_based,
    // x_final_0_based = img_xy_linear_absolute.x() + self.sip.f(u, v)
    let linear_x = img_xy_linear_absolute.x;
    let linear_y = img_xy_linear_absolute.y;
    let sip_offset_x = self.sip.f(u, v);
    let sip_offset_y = self.sip.g(u, v);

    eprintln!("[mapproj DEBUG] linear_x: {:.15}", linear_x);
    eprintln!("[mapproj DEBUG] sip_offset_x (A_u_v): {:.15}", sip_offset_x);
    eprintln!("[mapproj DEBUG] linear_y: {:.15}", linear_y);
    eprintln!("[mapproj DEBUG] sip_offset_y (B_u_v): {:.15}", sip_offset_y);

    let final_x = linear_x + sip_offset_x;
    let final_y = linear_y + sip_offset_y;
    eprintln!("[mapproj DEBUG] final_x (sum): {:.15}", final_x);
    eprintln!("[mapproj DEBUG] final_y (sum): {:.15}", final_y);
    Some(ImgXY::new(final_x, final_y))
  }
}
