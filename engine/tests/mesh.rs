use heatphase::thermal::material::{SILICA, SILICON};
use heatphase::thermal::mesh::{Mesh, Rect};

fn mesh() -> Mesh {
    Mesh::new(8.0, -2.0, 2.0, 0.5, 0.25, SILICA).unwrap()
}

#[test]
fn sizes_and_coordinates() {
    let m = mesh();
    assert_eq!((m.nx, m.ny, m.len()), (16, 16, 256));
    assert!(!m.is_empty());
    assert_eq!((m.x(0), m.y(0)), (-3.75, -1.875));
    assert_eq!(m.index(3, 2), 35);
}

#[test]
fn windows_must_be_whole_numbers_of_cells() {
    assert!(Mesh::new(8.1, -2.0, 2.0, 0.5, 0.25, SILICA).unwrap_err().contains("in x"));
    assert!(Mesh::new(8.0, -2.0, 2.1, 0.5, 0.25, SILICA).unwrap_err().contains("in y"));
    assert!(Mesh::new(8.0, -2.0, 2.0, 0.0, 0.25, SILICA).is_err());
    assert!(Mesh::new(0.5, -2.0, 2.0, 0.5, 0.25, SILICA).is_err());
}

#[test]
fn painting_goes_by_cell_centres() {
    let mut m = mesh();
    m.paint(&Rect { xmin: -1.0, xmax: 1.0, ymin: 0.0, ymax: 0.5 }, SILICON);
    let painted = m.cells.iter().filter(|c| **c == SILICON).count();
    assert_eq!(painted, 4 * 2);
    assert_eq!(m.cells[m.index(6, 8)], SILICON);
    assert_eq!(m.cells[m.index(5, 8)], SILICA);
    // a rectangle that contains no centre paints nothing
    m.paint(&Rect { xmin: 2.05, xmax: 2.2, ymin: -2.0, ymax: 2.0 }, SILICON);
    assert_eq!(m.cells.iter().filter(|c| **c == SILICON).count(), painted);
}

#[test]
fn coverage_adds_up_to_the_area_inside_the_window() {
    let m = mesh();
    let r = Rect { xmin: -0.7, xmax: 0.9, ymin: -0.1, ymax: 0.3 };
    let covered: f64 = m.coverage(&r).iter().map(|c| c.1).sum();
    assert!((covered - r.area()).abs() < 1e-12);
    // partly outside: only the inside counts
    let edge = Rect { xmin: 3.5, xmax: 5.0, ymin: 1.5, ymax: 3.0 };
    let inside: f64 = m.coverage(&edge).iter().map(|c| c.1).sum();
    assert!((inside - 0.25).abs() < 1e-12);
    assert!(m.coverage(&Rect { xmin: 10.0, xmax: 11.0, ymin: 0.0, ymax: 1.0 }).is_empty());
    assert_eq!(Rect { xmin: 1.0, xmax: 0.0, ymin: 0.0, ymax: 1.0 }.area(), 0.0);
}

#[test]
fn locating_points() {
    let m = mesh();
    assert_eq!(m.locate(-3.9, -1.9), Some(0));
    assert_eq!(m.locate(0.1, 0.1), Some(m.index(8, 8)));
    assert_eq!(m.locate(4.1, 0.0), None);
    assert_eq!(m.locate(0.0, -2.5), None);
}
