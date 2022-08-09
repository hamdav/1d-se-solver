use ndarray::{Array1, Array2, Axis, ArrayView1};
use ndarray::{concatenate, s, stack};
use ndarray_linalg::{Eigh, UPLO};
use std::f64::consts::PI;
use super::numerov::State;


pub fn find_bound_states(xbounds: (f64, f64),
                         support: (f64, f64), 
                         potential: &Vec<f64>) -> Vec<State> 
{
    /*
     * Find bound states using exact diagonalization
     *
     * xbounds is taken to be the places where the potential turns infinite.
     * support is where the potential is non-zero
     */


    // d is the number of basis states I will use.
    let d = 100;

    let hamiltonian = create_hamiltonian(xbounds, support, potential, d);

    // find the eigenstates
    let (eigvals, eigvecs) = hamiltonian.eigh(UPLO::Upper).unwrap();

    // Construct basis states (ALREADY DONE EARLIER...)
    let basis_states = calculate_basis_states(xbounds, potential.len(), d);

    // Construct the states, but only the ones with eigenvalue lower than 0
    eigvals.into_iter()
        .take_while(|&e| e < 0.)
        .zip(eigvecs.columns())
        .map(|(eigval, eigvec)| 
             State{
                 energy: eigval,
                 wf: eigvec.dot(&basis_states).to_vec(), // maybe use into vec instead for speed?
                 x_left: xbounds.0,
                 x_right: xbounds.1,
             })
        .collect()
}

fn calculate_basis_states(xbounds: (f64, f64), n: usize, d: usize) 
    -> Array2<f64>
{
    /*
     * Returns the basis states as the rows of a 2D array
     */
    let f = PI / (xbounds.1 - xbounds.0);
    let xs = Array1::linspace(xbounds.0, xbounds.1, n);
    // Array2::from((1..d)
    //     .map(|n| (xs.clone() * n as f64 * f).map(|t| t.sin()))
    //     .collect::<Vec<Array1<f64>>>())
    let mut rv = Array2::zeros((n, d));
    for mut row in rv.axis_iter_mut(Axis(0)) {
        // Perform calculations and assign to `row`; this is a trivial example:
        row.assign(&(&xs * n as f64 * f).map(|t| t.sin()));
    }
    rv
}



    

fn create_hamiltonian(xbounds: (f64, f64),
                      support: (f64, f64), 
                      potential: &Vec<f64>,
                      d: usize) -> Array2<f64>
{

    // Step one: construct the hamiltonian matrix.
    // This is done in the basis sin(n f (x-xbounds.0)), where f = pi / L
    // thus, ‹n|H|m› = 
    // \int sin(n f (x-x0)) 1/2 (mf)^2 (-sin(m f (x-x0))) 
    // + sin(n f (x-x0)) V(x) sin(m f (x-x0)) dx
    // = \delta_{nm} (mf)^2/2 + \int sin(n f (x-x0)) V(x) sin(m f (x-x0))
    //
    // Construct the full potential: add zeros to the beginning and end so that it
    // fills the xbounds
    let dx = (support.1 - support.0) / potential.len() as f64;
    let potential = concatenate(Axis(0), 
        &[Array1::zeros((((support.0 - xbounds.0) / dx).round() as usize,)).view(),
        Array1::from(potential.clone()).view(),
        Array1::zeros((((xbounds.1 - support.1) / dx).round() as usize,)).view(),])
        .unwrap();
    let xs = Array1::linspace(xbounds.0, xbounds.1, potential.len());

    let f = PI / (xbounds.1 - xbounds.0);
    let mut hamiltonian = Array2::from_diag(
        &Array1::range(1., 1. + d as f64, 1.)
        .map(|n| (n * f).powi(2) / 2.)
        );

    let basis_states = calculate_basis_states(xbounds, potential.len(), d);

    for n in 0..d {
        for m in n..d {
            hamiltonian[[n, m]] += (&basis_states.slice(s![n,..]) * &basis_states.slice(s![m,..]) * &potential).sum();
            // Maybe this part is actually unneccesary, see definition of eigh
            if m != n {
                hamiltonian[[m, n]] += hamiltonian[[n, m]];
            }
        }
    }
    hamiltonian
}
