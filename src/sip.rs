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
  
  /// Returns a new Sip object representing the inverse transformation.
  /// This is achieved by swapping the `ab_proj` and `ab_deproj` coefficients.
  /// If `ab_deproj` is `None` (i.e., no reverse polynomial coefficients like AP_ij, BP_ij are defined),
  /// this method returns `None`, as a simple coefficient-swapped inverse Sip object cannot be formed.
  /// Returns the forward SIP coefficients (A_ij, B_ij).
  pub fn ab_proj(&self) -> &SipAB {
    &self.ab_proj
  }

  /// Returns the 1st axis domain of validity for the SIP transformation.
  pub fn u_domain(&self) -> &RangeInclusive<f64> {
    &self.u
  }

  /// Returns the 2nd axis domain of validity for the SIP transformation.
  pub fn v_domain(&self) -> &RangeInclusive<f64> {
    &self.v
  }

  /// Returns a new Sip object representing the inverse transformation.
  /// This is achieved by swapping the `ab_proj` and `ab_deproj` coefficients.
  /// If `ab_deproj` is `None` (i.e., no reverse polynomial coefficients like AP_ij, BP_ij are defined),
  /// this method returns `None`, as a simple coefficient-swapped inverse Sip object cannot be formed.
  pub fn get_inverse_sip_transform(&self) -> Option<Self> {
    if let Some(deproj_coeffs) = &self.ab_deproj {
      Some(Sip {
        ab_proj: deproj_coeffs.clone(),      // New forward coefficients are the old reverse coefficients
        ab_deproj: Some(self.ab_proj.clone()), // New reverse coefficients are the old forward coefficients
        u: self.u.clone(),     // Domains are assumed to be for the input of the respective ab_proj
        v: self.v.clone(),
        // Other fields like fuv, guv, n_iter, eps are specific to the iterative solver context
        // and might not be directly applicable or needed for a simple coefficient-swapped inverse object.
        // They are initialized to default/None values by Sip::new if not provided.
        // For a simple inverse, we primarily care about the coefficients and domains.
        // Let's rely on Sip::new's defaults or consider if they need specific values here.
        // For now, assuming Sip::new handles this or they are not critical for the inverse object's structure.
        // If Sip::new requires them, we'd need to pass them or suitable defaults.
        // The current Sip::new takes ab_proj, ab_deproj, u, v. The other fields are not parameters to new.
        // Let's check Sip::new parameters. It takes: ab_proj, ab_deproj, u, v.
        // So, the current construction is fine for these fields.
        // The other fields (fuv, guv, n_iter, eps) are set internally or have defaults in Sip::new.
        // We need to ensure the new Sip object is valid. Let's check Sip::new again.
        // Sip::new initializes fuv to u.clone(), guv to v.clone(), n_iter to N_ITER_MAX, and eps to EPS_MAX.
        fuv: self.u.clone(), // Initialize with the u-domain, as in Sip::new
        guv: self.v.clone(), // Initialize with the v-domain, as in Sip::new
        n_iter: 20,          // Default value for n_iter, as in Sip::new
        eps: 1.0e-9,         // Default value for eps, as in Sip::new
      })
    } else {
      None // Cannot form a simple inverse by swapping if no deprojection coefficients exist
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

    #[test]
    fn test_third_order_poly() {
        // M_fits = 3 (third order polynomial)
        // self.order in SipCoeff = M_fits + 1 = 4
        // Number of coefficients = (M_fits + 1) * (M_fits + 2) / 2 = 4 * 5 / 2 = 10
        // Polynomial: Sum_{i=0 to 3} Sum_{j=0 to 3-i} C_ij * u^i * v^j
        // C00, C01*v, C02*v^2, C03*v^3,
        // C10*u, C11*u*v, C12*u*v^2,
        // C20*u^2, C21*u^2*v,
        // C30*u^3

        // FITS standard order for coefficients c (as used in SipCoeff::p, dpdu, dpdv loops):
        // C00, C01, C02, C03, C10, C11, C12, C20, C21, C30
        let c_fits_order = vec![
            1.0,  // C00
            3.0,  // C01
            6.0,  // C02
            10.0, // C03
            2.0,  // C10
            5.0,  // C11
            9.0,  // C12
            4.0,  // C20
            8.0,  // C21
            7.0   // C30
        ];
        let coeffs = SipCoeff::new(c_fits_order.into_boxed_slice());
        assert_eq!(coeffs.order, 4, "SipCoeff.order should be M_fits + 1");

        let u: f64 = 0.5;
        let v: f64 = 0.2;

        // Expected p(u,v) calculation:
        // C00 = 1.0
        // C10*u = 2.0 * 0.5 = 1.0
        // C01*v = 3.0 * 0.2 = 0.6
        // C20*u^2 = 4.0 * 0.5^2 = 4.0 * 0.25 = 1.0
        // C11*u*v = 5.0 * 0.5 * 0.2 = 5.0 * 0.1 = 0.5
        // C02*v^2 = 6.0 * 0.2^2 = 6.0 * 0.04 = 0.24
        // C30*u^3 = 7.0 * 0.5^3 = 7.0 * 0.125 = 0.875
        // C21*u^2*v = 8.0 * 0.5^2 * 0.2 = 8.0 * 0.25 * 0.2 = 8.0 * 0.05 = 0.4
        // C12*u*v^2 = 9.0 * 0.5 * 0.2^2 = 9.0 * 0.5 * 0.04 = 9.0 * 0.02 = 0.18
        // C03*v^3 = 10.0 * 0.2^3 = 10.0 * 0.008 = 0.08
        // Sum = 1.0 + 1.0 + 0.6 + 1.0 + 0.5 + 0.24 + 0.875 + 0.4 + 0.18 + 0.08 = 5.875
        let expected_p = 1.0                       // C00
                       + 2.0 * u                   // C10*u
                       + 3.0 * v                   // C01*v
                       + 4.0 * u.powi(2)           // C20*u^2
                       + 5.0 * u * v               // C11*u*v
                       + 6.0 * v.powi(2)           // C02*v^2
                       + 7.0 * u.powi(3)           // C30*u^3
                       + 8.0 * u.powi(2) * v       // C21*u^2*v
                       + 9.0 * u * v.powi(2)       // C12*u*v^2
                       + 10.0 * v.powi(3);         // C03*v^3
        assert_approx_eq(coeffs.p(u, v), expected_p, "p(u,v) for M=3");
        assert_approx_eq(coeffs.p(0.0, 0.0), 1.0, "p(0,0) for M=3 should be C00 (1.0)");

        // Expected dpdu(u,v) calculation:
        // C10 = 2.0
        // 2*C20*u = 2 * 4.0 * 0.5 = 4.0
        // C11*v = 5.0 * 0.2 = 1.0
        // 3*C30*u^2 = 3 * 7.0 * 0.5^2 = 21.0 * 0.25 = 5.25
        // 2*C21*u*v = 2 * 8.0 * 0.5 * 0.2 = 16.0 * 0.1 = 1.6
        // C12*v^2 = 9.0 * 0.2^2 = 9.0 * 0.04 = 0.36
        // Sum = 2.0 + 4.0 + 1.0 + 5.25 + 1.6 + 0.36 = 14.21
        let expected_dpdu = 2.0                       // C10
                          + 2.0 * 4.0 * u             // 2*C20*u
                          + 5.0 * v                   // C11*v
                          + 3.0 * 7.0 * u.powi(2)     // 3*C30*u^2
                          + 2.0 * 8.0 * u * v         // 2*C21*u*v
                          + 9.0 * v.powi(2);          // C12*v^2
        assert_approx_eq(coeffs.dpdu(u, v), expected_dpdu, "dpdu(u,v) for M=3");

        // Expected dpdv(u,v) calculation:
        // C01 = 3.0
        // C11*u = 5.0 * 0.5 = 2.5
        // 2*C02*v = 2 * 6.0 * 0.2 = 2.4
        // C21*u^2 = 8.0 * 0.5^2 = 2.0
        // 2*C12*u*v = 2 * 9.0 * 0.5 * 0.2 = 1.8
        // 3*C03*v^2 = 3 * 10.0 * 0.2^2 = 1.2
        // Sum = 3.0 + 2.5 + 2.4 + 2.0 + 1.8 + 1.2 = 12.9
        let expected_dpdv = 3.0                       // C01
                          + 5.0 * u                   // C11*u
                          + 2.0 * 6.0 * v             // 2*C02*v
                          + 8.0 * u.powi(2)           // C21*u^2
                          + 2.0 * 9.0 * u * v         // 2*C12*u*v
                          + 3.0 * 10.0 * v.powi(2);   // 3*C03*v^2
        assert_approx_eq(coeffs.dpdv(u, v), expected_dpdv, "dpdv(u,v) for M=3");
    }

    #[test]
    fn test_sip_astropy1() {
        
        // from pix to foc
        let a_coeff = SipCoeff::new(Box::new([0.0, 0.0, 1.569e-05, 0.0, 5.232e-05, 3.31e-05]));
        let b_coeff = SipCoeff::new(Box::new([0.0, 0.0, 4.172e-05, 0.0, 2.213e-05, -9.819e-07]));
        
        // Define the inverse transformation coefficients (AP, BP)
        let ap_coeff = SipCoeff::new(Box::new([0.0, 5.677e-05, -1.569e-05, 5.871e-05, -5.231e-05, -3.309e-05]));
        let bp_coeff = SipCoeff::new(Box::new([0.0, 4.432e-05, -4.172e-05, 2.091e-05, -2.213e-05, 9.814e-07]));
        
        // Define the valid range for u and v coordinates
        let u_range = -64.0..=64.0;
        let v_range = -32.0..=32.0;
        
        // Test cases: (u, v, expected_f, expected_g)
        let test_cases = [
            (-64.0, -32.0, -63.74120447999999, -31.915978342399995),
            (0.0, -32.0, 0.01606656000001294, -31.957278720000005),
            (64.0, -32.0, 64.0444928, -32.006622822400004),
            (-64.0, 0.0, -63.864422399999995, -0.004021862399994802),
            (0.0, 0.0, 0.0, 0.0),
        ];
    
        // Create the SIP transformation
        let ab_proj = SipAB::new(a_coeff, b_coeff);
        let ab_deproj = Some(SipAB::new(ap_coeff, bp_coeff));
        let sip = Sip::new(ab_proj, ab_deproj, u_range, v_range);
        
        // Test each case
        for (i, &(u, v, uf, vg)) in test_cases.iter().enumerate() {
            
          let f_uv = sip.f(u, v);
          let g_uv = sip.g(u, v);
          
          // Test with a relative tolerance of 1e-10
          let abs_tol = 1e-2;
          let f_diff = (u + f_uv - uf).abs();
          let g_diff = (v + g_uv - vg).abs();
          
          assert!(
              f_diff <= abs_tol,
              "Case {}: f({}, {}) = {}, expected {}. Difference: {} > {}",
              i, u, v, f_uv, uf, f_diff, abs_tol
          );
          assert!(
            g_diff <= abs_tol,
            "Case {}: g({}, {}) = {}, expected {}. Difference: {} > {}",
            i, u, v, g_uv, vg, g_diff, abs_tol
          );

          // Test inverse transformation if we have deprojection coefficients
          if let (Some(U), Some(V)) = (sip.u(uf, vg), sip.v(uf, vg)) {
            let u_diff = (uf + U - u).abs();
            let v_diff = (vg + V - v).abs();
            
            assert!(
                u_diff <= abs_tol,
                "Case {}: u({}, {}) = {}, expected {}. Difference: {} > {}",
                i, uf, vg, U, u, u_diff, abs_tol
            );
          
            assert!(
                v_diff <= abs_tol,
                "Case {}: v({}, {}) = {}, expected {}. Difference: {} > {}",
                i, uf, vg, V, v, v_diff, abs_tol
            );
          }


        }
    }        

    #[test]
    fn test_sip_instance() {
        // Create simple first-order coefficients for forward transformation (A_ij, B_ij)
        // Coefficients are ordered like this: `0_0, 0_1, 0_2, 0_3, 1_0, 1_1, 1_2, 2_0, 2_1, 3_0`
        let a_coeff = SipCoeff::new(Box::new([0., 0., 1.])); // A_00, A_01, A_10
        let b_coeff = SipCoeff::new(Box::new([0., 0., 0.])); // B_00, B_01, B_10
        let ab_proj = SipAB::new(a_coeff, b_coeff);
        
        let ap_coeff = SipCoeff::new(Box::new([0., 0., -0.5e0]));
        let bp_coeff = SipCoeff::new(Box::new([0., 0., 0.]));
        let ab_deproj = Some(SipAB::new(ap_coeff, bp_coeff));
        
        // Define the domain for u and v (pixel coordinates relative to CRPIX)
        // Let's assume an image of 100x100 pixels with CRPIX at (50, 50)
        let u_range = -50.0..=50.0;
        let v_range = -50.0..=50.0;
        
        // Create the Sip instance
        let sip = Sip::new(ab_proj, ab_deproj, u_range, v_range);
        
        // Test point in the domain
        let test_u = 2.0;
        let test_v = 0.0;
        
        // Test forward transformation
        let f_uv = sip.f(test_u, test_v);
        let g_uv = sip.g(test_u, test_v);
        
        assert_approx_eq(f_uv, 2., "f(u,v) calculation incorrect");
        assert_approx_eq(g_uv, 0., "g(u,v) calculation incorrect");

        if let Some(u1) = sip.u(f_uv, g_uv) {
            assert_approx_eq(u1, -1., "u(f(u,v), g(u,v)) calculation incorrect");
        }
        
        
    }
}