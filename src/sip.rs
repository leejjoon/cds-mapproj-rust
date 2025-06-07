//! Implementation of the SIP standard.

use std::ops::RangeInclusive;
use crate::{CustomFloat, ImgXY};

/// SIP Polynomial coefficient.
/// In the polynomial, coefficient must be ordered like this:
/// * `0_0, 0_1, 0_2, 0_3, 1_0, 1_1, 1_2, 2_0, 2_1, 3_0`
/// * in which `p_q` correspond to the polynomial part `coeff_p_q * u^p * v^q` 
/// * Given an order `n`, the size of the array must be `n(n+1)/2`.
#[derive(Clone)]
pub struct SipCoeff {
  /// Computed order of the polynomial
  order: u16,
  /// Polynomials coefficient matrix
  c: Box<[f64]>,
}

impl SipCoeff {
  /// # Param
  /// * `c`: array polynomial coefficients of size `n(n+1)/2`
  pub fn new(c: Box<[f64]>) -> Self {
    let order = Self::order_from_n_coeff(c.len());
    debug_assert_eq!(order * (order + 1), (c.len() as u16) << 1);
    Self { order, c }
  }
  
  /// Returns the order of a bivariate polynomial from its number of coefficient.
  /// Sum of k for k = 1 to n equals n(n+1)/2 = l
  /// Thus n = (sqrt(8*l + 1) - 1) / 2
  /// # Param 
  /// * `n_coeff`: number of coefficient of the bivariate polynomial
  fn order_from_n_coeff(n_coeff: usize) -> u16 {
    ((((n_coeff << 3) + 1) as f64).sqrt()  as u16 - 1) >> 1
  }
  
  /// Returns the value of the polynomial, evaluated in `(u, v)`.
  pub fn p(&self, u: f64, v: f64) -> f64 {
    if self.c.is_empty() { // Handle case of no coefficients
        return 0.0;
    }
    let mut sum_val = 0.0;
    let mut k = 0; // Index for the coefficient array self.c
    // self.order stores M_fits + 1, where M_fits is the FITS polynomial degree.
    let m_fits = self.order - 1;
    // Polynomial: Sum_{i=0 to M_fits} Sum_{j=0 to M_fits-i} C_ij * u^i * v^j
    for i_power in 0..=m_fits { // i from 0 to M_fits (inclusive)
        let u_pow_i = u.powi(i_power as i32);
        for j_power in 0..=(m_fits - i_power) { // j from 0 to M_fits-i (inclusive)
            // Given the loop structure and m_fits derivation, k should not exceed c.len()
            // before the loop finishes all terms.
            let v_pow_j = v.powi(j_power as i32);
            sum_val += self.c[k] * u_pow_i * v_pow_j;
            k += 1;
        }
    }
    debug_assert_eq!(k, self.c.len(), "Mismatch in coefficient usage for p: expected {}, got {}", self.c.len(), k);
    sum_val
  }

  /// Returns the value of the `dp/du`, evaluated in `(u, v)`.
  pub fn dpdu(&self, u: f64, v: f64) -> f64 {
    if self.c.is_empty() { return 0.0; }
    let mut sum_val = 0.0;
    let mut k = 0;
    let m_fits = self.order - 1;
    // Derivative: Sum C_ij * i * u^(i-1) * v^j
    for i_power in 0..=m_fits {
        for j_power in 0..=(m_fits - i_power) {
            if i_power > 0 {
                let u_pow_i_minus_1 = u.powi(i_power as i32 - 1);
                let v_pow_j = v.powi(j_power as i32);
                sum_val += self.c[k] * (i_power as f64) * u_pow_i_minus_1 * v_pow_j;
            }
            k += 1;
        }
    }
    debug_assert_eq!(k, self.c.len(), "Mismatch in coefficient usage for dpdu: expected {}, got {}", self.c.len(), k);
    sum_val
  }

  /// Returns the value of the `dp/dv`, evaluated in `(u, v)`.
  pub fn dpdv(&self, u: f64, v: f64) -> f64 {
    if self.c.is_empty() { return 0.0; }
    let mut sum_val = 0.0;
    let mut k = 0;
    let m_fits = self.order - 1;
    // Derivative: Sum C_ij * u^i * j * v^(j-1)
    for i_power in 0..=m_fits {
        for j_power in 0..=(m_fits - i_power) {
            if j_power > 0 {
                let u_pow_i = u.powi(i_power as i32);
                let v_pow_j_minus_1 = v.powi(j_power as i32 - 1);
                sum_val += self.c[k] * (j_power as f64) * u_pow_i * v_pow_j_minus_1;
            }
            k += 1;
        }
    }
    debug_assert_eq!(k, self.c.len(), "Mismatch in coefficient usage for dpdv: expected {}, got {}", self.c.len(), k);
    sum_val
  }
  
}


/// SIP (un)projection coefficients for 1st and 2nd axis
#[derive(Clone)]
pub struct SipAB {
  /// Polynomials coefficient matrix on the 1st axis.
  a: SipCoeff,
  /// Polynomials coefficient matrix on the2ndt axis.
  b: SipCoeff
}

impl SipAB {
  /// # Params
  /// * `a`: 1st axis SIP coefficients
  /// * `b`: 2nd axis SIP coefficients
  pub fn new(a: SipCoeff, b: SipCoeff) -> Self {
    Self { a, b }
  }
}

/// For the SIP convention, see
/// "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
/// in the proceedings of ADASS XIV (2005).
#[derive(Clone)]
pub struct Sip {
  /// Projection coefficient.
  ab_proj: SipAB,
  /// Unprojection coefficients (if any).
  ab_deproj: Option<SipAB>, // Unprojection coefficients are optional :o/
  /// Approximatve bounds of the 1st axis domain of validity.
  /// * `start`: `-(CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS1 - CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  u: RangeInclusive<f64>,
  /// Approximatve bounds of the 2nd axis domain of validity.
  /// * `start`: `-(CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS2 - CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  v: RangeInclusive<f64>,
  fuv: RangeInclusive<f64>,
  guv: RangeInclusive<f64>,
  /// Number of iteration of the mutli-variate Newton-Raphson method (if no unproj polynomial).
  n_iter: u8, // = 20;
  /// Precision used in the mutli-variate Newton-Raphson method (if no unproj polynomial).
  eps: f64,   // = 1e-9;
}

impl Sip {
  
  /// Implements the SIP convention with the given polynomial coefficients.
  /// # Params
  /// * `ab_proj`: SIP coefficients for the projection on the 1st and 2nd axis 
  /// * `ab_deproj`: SIP coefficients for the deprojection on the 1st and 2nd axis (if any)
  /// * `u`: 1st axis domain of validity, e.g. `[-CRPIX1..NAXIS1 - CRPIX1]` 
  /// * `v`: 2nd axis domain of validity, e.g. `[-CRPIX2..NAXIS2 - CRPIX2]`
  pub fn new(
    ab_proj: SipAB, 
    ab_deproj: Option<SipAB>, 
    u: RangeInclusive<f64>, 
    v: RangeInclusive<f64>,
  ) -> Self {
    let t = ab_proj.a.p(*u.start(), *v.start()).min(ab_proj.a.p(*u.start(), *v.end()));
    let fuv_min = ab_proj.a.p(*u.start(), 0.0).min(t);
    let t = ab_proj.a.p(*u.end(), *v.start()).max(ab_proj.a.p(*u.end(), *v.end()));
    let fuv_max = ab_proj.a.p(*u.end(), 0.0).max(t);
    
    let t =  ab_proj.b.p(*u.start(), *v.start()).min( ab_proj.b.p(*u.end(), *v.start()));
    let guv_min =  ab_proj.b.p(0.0, *v.start()).min(t);
    let t =  ab_proj.b.p(*u.start(), *v.end()).max( ab_proj.b.p(*u.end(), *v.end()));
    let guv_max =  ab_proj.b.p(0.0, *v.end()).max(t);

    Self { 
      ab_proj, 
      ab_deproj, 
      u, v,
      fuv: fuv_min..=fuv_max,
      guv: guv_min..=guv_max,
      n_iter: 20,
      eps: 1.0e-9
    }
  }
  
  pub fn has_polynomial_deproj(&self) -> bool {
    self.ab_deproj.is_some()
  }
  
  pub fn f(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.p(u, v)
  }

  pub fn g(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.p(u, v)
  }

  pub fn dfdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdu(u, v)
  }

  pub fn dfdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdv(u, v)
  }

  pub fn dgdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdu(u, v)
  }

  pub fn dgdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdv(u, v)
  }

  /// Evaluates the AP polynomial (inverse distortion for the first axis).
  /// Arguments u, v are 1-indexed pixel offsets: (pixel_coord - CRPIX_header).
  pub fn ap(&self, u: f64, v: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|coeffs| coeffs.a.p(u, v))
  }

  /// Evaluates the BP polynomial (inverse distortion for the second axis).
  /// Arguments u, v are 1-indexed pixel offsets: (pixel_coord - CRPIX_header).
  pub fn bp(&self, u: f64, v: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|coeffs| coeffs.b.p(u, v))
  }

  pub fn u(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.a.p(fuv, guv))
  }

  pub fn v(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.b.p(fuv, guv))
  }

  pub fn inverse(&self, fuv: f64, guv: f64) -> Option<ImgXY> { // uv
    if self.has_polynomial_deproj() {
      let u = self.u(fuv, guv).unwrap();
      let v = self.v(fuv, guv).unwrap();
      Some(ImgXY::new(u, v))
    } else {
      // Make a grid and a 2-d tree to find the starting point (then multi-variate Newton) 
      None
    }
  }
  
  /// Mutli-variate Newton-Raphson:
  /// f1(x1, ..., xn) = 0es006500
  /// ...
  /// fn(x1, ..., xn) = 0
  ///  
  /// x = (x1, ..., xn)
  /// f = (f1(x), ... fn(x))
  ///  
  /// x = x - J^-1 f
  ///
  ///  With J = df1/dx1 ... df1/dxn
  /// ...     ... ...
  /// dfn/dx1 ... dfn/dxn
  ///
  /// 2d case: M = ab => M^-1 = 1/(ad-bc)  d -b 
  /// cd                     -c  a
  ///
  pub fn bivariate_newton(&self, fuv :f64, guv: f64) -> Option<ImgXY> {
    // Check input values are in the domain of validity
    if self.fuv.contains(&fuv) && self.guv.contains(&guv) {
      // Initial guess
      let mut u = fuv;
      let mut v = guv;
      // Initial values
      let mut f = self.f(u, v) - fuv;
      let mut g = self.g(u, v) - guv;
      // Bivariate Newton's method
      let eps2 = self.eps.pow2();
      let mut norm2 = f.pow2() + g.pow2();
      let mut i = 0;
      while i < self.n_iter && norm2 < eps2 {
        let a = self.dfdu(u, v);
        let b = self.dfdv(u, v);
        let c = self.dgdu(u, v);
        let d = self.dgdv(u, v);
        let det = 1.0 / (a * d - b * c);
        u -= det * (f * d - g * b);
        v -= det * (g * a - f * c);
        f = self.f(u, v) - fuv;
        g = self.g(u, v) - guv;
        norm2 = f.pow2() + g.pow2();
        i += 1;
      }
      // Check that the result is in the domain of validity
      if self.u.contains(&u) && self.v.contains(&v) {
        // All good :)
        Some(ImgXY::new(u, v))
      } else {
        // TODO: look for a different initial guess by making a grid of 300x300 an put results
        // in a kd-tree. Then look at the nearest neighbour: it is the initial guess.
        // And redo Newton.
        None
      }
    } else {
      None
    }
  }
}

#[cfg(test)]
mod sip_coeff_tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn assert_approx_eq(a: f64, b: f64, msg: &str) {
        assert!((a - b).abs() < EPS, "{}\n  left: {}\n right: {}", msg, a, b);
    }

    #[test]
    fn test_empty_coeffs() {
        let coeffs = SipCoeff::new(Box::new([]));
        assert_eq!(coeffs.order, 0, "Order for empty coeffs should be 0");
        assert_approx_eq(coeffs.p(1.0, 1.0), 0.0, "p(u,v) with empty coeffs");
        assert_approx_eq(coeffs.dpdu(1.0, 1.0), 0.0, "dpdu(u,v) with empty coeffs");
        assert_approx_eq(coeffs.dpdv(1.0, 1.0), 0.0, "dpdv(u,v) with empty coeffs");
    }

    #[test]
    fn test_zero_order_poly() {
        // M_fits = 0. f(u,v) = C00
        // self.order = M_fits + 1 = 1
        let c = vec![10.0]; // C00 = 10.0
        let coeffs = SipCoeff::new(c.into_boxed_slice());
        assert_eq!(coeffs.order, 1);

        let u = 2.0;
        let v = 3.0;

        // p(u,v) = C00
        assert_approx_eq(coeffs.p(u, v), 10.0, "p(u,v) for M=0");
        assert_approx_eq(coeffs.p(0.0, 0.0), 10.0, "p(0,0) for M=0");
        // dpdu(u,v) = 0
        assert_approx_eq(coeffs.dpdu(u, v), 0.0, "dpdu(u,v) for M=0");
        // dpdv(u,v) = 0
        assert_approx_eq(coeffs.dpdv(u, v), 0.0, "dpdv(u,v) for M=0");
    }

    #[test]
    fn test_first_order_poly() {
        // M_fits = 1. f(u,v) = C00 + C10*u + C01*v
        // self.order = M_fits + 1 = 2
        // Conceptual coefficient values: C00=1.0, C10=2.0, C01=3.0
        // FITS standard order for c: C00, C01, C10
        let c = vec![1.0, 3.0, 2.0];
        let coeffs = SipCoeff::new(c.into_boxed_slice());
        assert_eq!(coeffs.order, 2);

        let u = 0.5;
        let v = -0.2;

        // p(u,v) = 1.0 + 2.0*0.5 + 3.0*(-0.2) = 1.0 + 1.0 - 0.6 = 1.4
        let expected_p = 1.0 + 2.0 * u + 3.0 * v;
        assert_approx_eq(coeffs.p(u, v), expected_p, "p(u,v) for M=1");
        assert_approx_eq(coeffs.p(0.0, 0.0), 1.0, "p(0,0) for M=1 should be C00");

        // dpdu(u,v) = C10 = 2.0
        let expected_dpdu = 2.0;
        assert_approx_eq(coeffs.dpdu(u, v), expected_dpdu, "dpdu(u,v) for M=1");

        // dpdv(u,v) = C01 = 3.0
        let expected_dpdv = 3.0;
        assert_approx_eq(coeffs.dpdv(u, v), expected_dpdv, "dpdv(u,v) for M=1");
    }

    #[test]
    fn test_second_order_poly() {
        // M_fits = 2. f(u,v) = C00 + C10*u + C01*v + C20*u^2 + C11*u*v + C02*v^2
        // self.order = M_fits + 1 = 3
        // Conceptual coefficient values: C00=1.0, C10=2.0, C01=3.0, C20=4.0, C11=5.0, C02=6.0
        // FITS standard order for c: C00, C01, C02, C10, C11, C20
        let c = vec![1.0, 3.0, 6.0, 2.0, 5.0, 4.0];
        let coeffs = SipCoeff::new(c.into_boxed_slice());
        assert_eq!(coeffs.order, 3);

        let u: f64 = 0.5;
        let v: f64 = 0.2;

        // p(u,v) = 1.0 + 2.0*u + 3.0*v + 4.0*u^2 + 5.0*u*v + 6.0*v^2
        // p(0.5, 0.2) = 1.0 + 2*0.5 + 3*0.2 + 4*0.5^2 + 5*0.5*0.2 + 6*0.2^2
        //             = 1.0 + 1.0   + 0.6   + 4*0.25  + 5*0.1   + 6*0.04
        //             = 1.0 + 1.0   + 0.6   + 1.0     + 0.5     + 0.24
        //             = 4.34
        let expected_p = 1.0 + 2.0*u + 3.0*v + 4.0*u.powi(2) + 5.0*u*v + 6.0*v.powi(2);
        assert_approx_eq(coeffs.p(u, v), expected_p, "p(u,v) for M=2");
        assert_approx_eq(coeffs.p(0.0, 0.0), 1.0, "p(0,0) for M=2 should be C00");

        // dpdu(u,v) = C10 + 2*C20*u + C11*v
        // dpdu(0.5, 0.2) = 2.0 + 2*4.0*0.5 + 5.0*0.2
        //                = 2.0 + 4.0       + 1.0
        //                = 7.0
        let expected_dpdu = 2.0 + 2.0*4.0*u + 5.0*v;
        assert_approx_eq(coeffs.dpdu(u, v), expected_dpdu, "dpdu(u,v) for M=2");

        // dpdv(u,v) = C01 + C11*u + 2*C02*v
        // dpdv(0.5, 0.2) = 3.0 + 5.0*0.5 + 2*6.0*0.2
        //                = 3.0 + 2.5       + 2.4
        //                = 7.9
        let expected_dpdv = 3.0 + 5.0*u + 2.0*6.0*v;
        assert_approx_eq(coeffs.dpdv(u, v), expected_dpdv, "dpdv(u,v) for M=2");
    }
}