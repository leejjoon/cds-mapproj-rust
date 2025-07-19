use mapproj::{
  zenithal::sin::Sin,
  CanonicalProjection, ProjXY, XYZ,
};

#[test]
fn test_sin_projection() {
  let sin = Sin::new();

  // Test projection
  let xyz = XYZ::new(1.0, 0.0, 0.0);
  let xy = sin.proj(&xyz).unwrap();
  assert!((xy.x() - 0.0).abs() < 1e-14);
  assert!((xy.y() - 0.0).abs() < 1e-14);

  let xyz = XYZ::new_renorming_if_necessary(0.5, 0.5, 0.5);
  let xy = sin.proj(&xyz).unwrap();
  let expected = 1.0 / 3.0_f64.sqrt();
  assert!((xy.x() - expected).abs() < 1e-14);
  assert!((xy.y() - expected).abs() < 1e-14);

  // Test projection of a point on the back hemisphere
  let xyz = XYZ::new_renorming_if_necessary(-0.5, 0.5, 0.5);
  assert!(sin.proj(&xyz).is_none());
}

#[test]
fn test_sin_unprojection() {
  let sin = Sin::new();

  // Test unprojection
  let xy = ProjXY::new(0.0, 0.0);
  let xyz = sin.unproj(&xy).unwrap();
  assert!((xyz.x() - 1.0).abs() < 1e-14);
  assert!((xyz.y() - 0.0).abs() < 1e-14);
  assert!((xyz.z() - 0.0).abs() < 1e-14);

  let xy = ProjXY::new(0.5, 0.5);
  let xyz = sin.unproj(&xy).unwrap();
  assert!((xyz.x() - (0.5_f64).sqrt()).abs() < 1e-14);
  assert!((xyz.y() - 0.5).abs() < 1e-14);
  assert!((xyz.z() - 0.5).abs() < 1e-14);

  // Test unprojection of a point outside the projection plane
  let xy = ProjXY::new(1.1, 0.0);
  assert!(sin.unproj(&xy).is_none());
}

#[test]
fn test_sin_consistency() {
    let sin = Sin::new();
    // Test projection and unprojection consistency
    let xyz_original = XYZ::new_renorming_if_necessary(0.5, 0.5, 0.5);
    let xy = sin.proj(&xyz_original).unwrap();
    let xyz_new = sin.unproj(&xy).unwrap();
    assert!((xyz_original.x() - xyz_new.x()).abs() < 1e-14);
    assert!((xyz_original.y() - xyz_new.y()).abs() < 1e-14);
    assert!((xyz_original.z() - xyz_new.z()).abs() < 1e-14);
}
