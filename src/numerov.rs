use super::Float;


#[derive(Debug)]
pub struct State{
    pub energy: Float,
    pub wf: Vec<Float>,
    pub x_left: Float,
    pub x_right: Float,
}

pub fn find_bound_states(xbounds: (Float, Float),
                         support: (Float, Float), 
                         potential: &Vec<Float>) -> Vec<State> 
{
    /*
     * Find the bound states using the numerov algorithm
     *
     * Assumes that the values of the potential are sampled on evenly spaced
     * x-values on the support (endpoints included).
     *
     * The wavefunctions of the bound state are sampled on evenly spaced
     * x-values on the xbounds (again, endpoints included)
     */

    // Find the minimum energy of the potential, no bound state can have
    // and energy smaller than this
    let E_min = *potential.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    // The largest energy a bound state can have is technically infinitely
    // close to 0. but since that runs into problems with normalization
    // I cap it at E_min / 1000. No one will be able to tell the difference
    // anyway
    let energy_interval = (E_min, E_min/1000.);

    // Find the wavefunctions at the energy interval endpoints, for use
    // as the initial values of the bisection algorithm
    let (psi_lo, _) = bidirectional_shooting(energy_interval.0, xbounds, 
                                             support, potential);
    let (psi_hi, _) = bidirectional_shooting(energy_interval.1, xbounds, 
                                             support, potential);

    // Return value. Will contain the energy eigen-wavefunctions
    let mut rv = Vec::new();

    // We know that there is exactly one energy eigenstate with
    // any given number of nodes up to some limit, and that the
    // more nodes the state has, the higher energy it has.
    // 
    // So our algorithm first finds an interval enclosing the boundry
    // between n and n+1 nodes, and then uses the normal bisection
    // method to find the energy where the log_diff is 0 in the interval
    // with n nodes. The lower bound for the 0 node interval is E_min
    // and the upper bound for the last interval is E_min / 1000.
    let mut lower_end = E_min;
    for n_nodes in 0.. {

        // Find interval containing the boundry between n and n+1 nodes
        let interval = find_interval_enclosing_upper_boundry(
            energy_interval,
            (count_nodes(&psi_lo), count_nodes(&psi_hi)),
            n_nodes, xbounds, support, potential);

        // If there was one, do the bisection
        if let Some(interval) = interval {

            let (psi_lo, log_diff_lo) = bidirectional_shooting(
                lower_end, xbounds, support, potential);
            let (psi_hi, log_diff_hi) = bidirectional_shooting(
                interval.0, xbounds, support, potential);

            println!("");
            println!("E: {:?}", (lower_end, interval.0));
            println!("log_diffs: {:?}", (log_diff_lo, log_diff_hi));
            println!("# of nodes: {:?}", (count_nodes(&psi_lo), count_nodes(&psi_hi)));
            let E = simple_logdiff_bisect((lower_end, interval.0),
                (log_diff_lo, log_diff_hi),
                xbounds, support, potential);

            // If the bisection was successful, add the wavefunction.
            // Otherwise, print an error message
            if let Some(E) = E {
                let (psi, _) = bidirectional_shooting(E, xbounds, support, potential);
                rv.push(psi);
            } else {
                println!("Failed logdiff bisection");
                // So the bisection failed. This probably means that the
                // eigen energy is closer to the node boundry than the size
                // of the interval. If so, one of the boundaries is
                // probably close enough to the true wf, the only question
                // is which. To determine this, we calculate the overlap
                // with the previously added wf and take the one that is
                // most orthogonal.
                let dx = (xbounds.1 - xbounds.0) / psi_lo.wf.len() as Float;
                let overlaps = (
                    psi_lo.wf.iter().zip(rv[rv.len()-1].wf.iter())
                        .map(|(&pl, &p)| pl * p)
                        .sum::<Float>() * dx,
                    psi_hi.wf.iter().zip(rv[rv.len()-1].wf.iter())
                        .map(|(&ph, &p)| ph * p)
                        .sum::<Float>() * dx
                    );

                println!("overlaps: {:?}", overlaps);
                if overlaps.0 < overlaps.1 {
                    rv.push(psi_lo);
                } else {
                    rv.push(psi_hi);
                }
            }

            // Prepare for next loop by saving the lower bound
            lower_end = interval.1;
        } else {
            // This will be the last wavefunction (if there is one)
            // Do one loop of the same as above with the upper end being 
            // E_min / 1000.
            let (psi_lo, log_diff_lo) = bidirectional_shooting(
                lower_end, xbounds, support, potential);
            let (psi_hi, log_diff_hi) = bidirectional_shooting(
                energy_interval.1, xbounds, support, potential);

            println!("");
            println!("E: {:?}", (lower_end, energy_interval.1));
            println!("log_diffs: {:?}", (log_diff_lo, log_diff_hi));
            println!("# of nodes: {:?}", (count_nodes(&psi_lo), count_nodes(&psi_hi)));
            let E = simple_logdiff_bisect((lower_end, energy_interval.1),
                (log_diff_lo, log_diff_hi),
                xbounds, support, potential);

            // No error message needs to be printed now because we don't
            // expect that this will succeed always.
            if let Some(E) = E {
                let (psi, _) = bidirectional_shooting(E, xbounds, support, potential);
                rv.push(psi);
            }
            break;
        }
    }
    rv
} 

fn count_nodes(psi: &State) -> usize {
    psi.wf.iter()
        .zip(psi.wf.iter().skip(1))
        .map(|(&a, &b)| if a * b < 0. {1} else {0})
        .sum::<usize>()
}

fn simple_logdiff_bisect(mut energy_interval: (Float, Float),
                         mut end_log_diffs: (Float, Float),
                         xbounds: (Float, Float),
                         support: (Float, Float),
                         potential: &Vec<Float>)
    -> Option<Float>
{
    /*
     * Returns Some(energy) if the different ends of the energy interval
     * has the same number of nodes and different signs for the log diff
     */


    // If the log diffs are both supersmall, return the middle
    // if end_log_diffs.0.abs() < 1e-3 {
    //     return Some(energy_interval.0)
    // } else if end_log_diffs.1.abs() < 1e-3 {
    //     return Some(energy_interval.1)
    // }

    // Check that interval is valid
    if end_log_diffs.0 * end_log_diffs.1 > 0. {
        return None
    }

    // Base case: interval is small = return middle of interval.
    let mid = (energy_interval.0 + energy_interval.1) / 2.;
    if ((energy_interval.1 - energy_interval.0) / energy_interval.0).abs() < 1e-13 {
        return Some(mid)
    }

    // Check mid
    let (_, log_diff) = bidirectional_shooting(mid, xbounds, support, potential);

    // Change appropriate bound
    if log_diff * end_log_diffs.0 <= 0. {
        end_log_diffs.1 = log_diff;
        energy_interval.1 = mid;
    } else {
        end_log_diffs.0 = log_diff;
        energy_interval.0 = mid;
    }
    // Recurse
    return simple_logdiff_bisect(energy_interval, end_log_diffs, xbounds, support, potential);
}

fn find_interval_enclosing_upper_boundry(mut energy_interval: (Float, Float),
                                   mut end_nodes: (usize, usize),
                                   n_nodes: usize,
                                   xbounds: (Float, Float),
                                   support: (Float, Float),
                                   potential: &Vec<Float>)
    -> Option<(Float, Float)>
{
    /*
     * Finds a small interval enclosing the upper endpoint of the energy 
     * interval where the wf has n_nodes nodes, if there is one.
     */

    // Base case: interval is small = return lower end of interval.
    if ((energy_interval.1 - energy_interval.0) / energy_interval.0).abs() < 1e-13 {
        if end_nodes.0 == n_nodes && end_nodes.1 == n_nodes+1 {
            return Some(energy_interval)
        } else {
            return None
        }
    }
    // Check that the interval is even valid
    if end_nodes.1 <= n_nodes {
        return None
    }

    // Do the bisection
    let mid = (energy_interval.0 + energy_interval.1) / 2.;
    let (psi, log_diff) = bidirectional_shooting(mid, xbounds, support, potential);
    let mid_n_nodes = count_nodes(&psi);
    // println!("");
    // println!("E: {:?}", energy_interval);
    // println!("# of nodes: {:?}", end_nodes);
    // println!("mid: {}, {}, {}", mid, mid_n_nodes, log_diff);

    // Change one of the bounds
    if mid_n_nodes <= n_nodes {
        end_nodes.0 = mid_n_nodes;
        energy_interval.0 = mid;
    } else {
        end_nodes.1 = mid_n_nodes;
        energy_interval.1 = mid;
    }

    // Recurse
    return find_interval_enclosing_upper_boundry(energy_interval, end_nodes, n_nodes, xbounds, support, potential)
}

fn min_index<I>(iter: I) -> usize
    where I: IntoIterator<Item = Float>
{
    /*
     * finds the index of a minimum of an IntoIterator over Float elements
     */
    let mut min = Float::INFINITY;
    let mut mindex = 0;

    for (i, v) in iter.into_iter().enumerate() {
        if v < min {
            min = v;
            mindex = i;
        }
    }
    mindex
}
                
fn numerov(E: Float, psi_0: Float, psi_1: Float, start_ind: usize, end_ind: usize, dx: Float,
           potential: &Vec<Float>) -> Vec<Float>
{
    /*
     * Performs numerov integration from index start_ind of the potential
     * to index end_ind (inclusive) given the values of the function
     * at the index before start_ind (psi_0) and at start_ind (psi_1).
     * Note that the index before start_ind may be start_ind+1 if end_ind is
     * smaller than start_ind, as we then go backwards.
     *
     * The function values returned are in the order of integration,
     * so if start_ind = 10 and end_ind = 5 the returned vector would be
     * [f(11), f(10), f(9), f(8), f(7), f(6), f(5)]
     *
     * Also note that you cannot start at 0, or potential.len()-1.
     */
    let mut psis = Vec::with_capacity(
        (end_ind as i32 - start_ind as i32).abs() as usize);
    psis.push(psi_0);
    psis.push(psi_1);

    // Create a list of indices to iterate over
    // would be start_ind..=end_ind except that fails if start_ind is
    // larger than end_ind.
    let inds: Vec<usize> =
        if start_ind < end_ind {
            (start_ind..=end_ind).collect()
        } else {
            (end_ind..=start_ind).rev().collect()
        };

    for i in inds {
        let a = 1. + dx.powi(2)/12. * 2.*(E - potential[i+1]);
        let b = 1. - 5.*dx.powi(2)/12. * 2.*(E - potential[i]);
        let c = 1. + dx.powi(2)/12. * 2.*(E - potential[i-1]);

        psis.push((2.*b*psis[psis.len()-1] - c*psis[psis.len() - 2]) / a);
    }

    psis
}
    


fn bidirectional_shooting(E: Float, 
                          xbounds: (Float, Float), 
                          support: (Float, Float), 
                          potential: &Vec<Float>)
    -> (State, Float)
{
    /*
     * runs numerov from left and from right to b with energy E
     * Potential is assumed to be sampled on evenly spaced points with
     * the first one being support.0 and the last support.1
     * The returned wavefunction will be sampled with the same dx as 
     * potential but on the interval defined by xbounds.
     *
     * TODO WHATIS B?
     *
     *
     * Returns the concatinated wavefunction along with
     */

    assert!(E < 0.);

    let dx = (support.1 - support.0) / (potential.len() - 1) as Float;

    //let mut i = min_index(potential.iter().map(|v| (v-E).abs()));
    let mut i = potential.iter()
        .skip_while(|&&v| v > E)
        .skip_while(|&&v| v < E)
        .count();
    i = potential.len() - i;
    if i <= 1 { i = 2; }
    else if i >= potential.len() - 2 {i = potential.len()-3;}

    // if V[i] - E is too large, do something... does it actually matter?
    // I don't think so
    let psi_l = numerov(E, (-(-2.*E).sqrt() * dx).exp(), 1., 1, i, dx, &potential);
    let psi_r = numerov(E, (-(-2.*E).sqrt() * dx).exp(), 1., potential.len()-2, i, dx, &potential);

    // Find the norms of the right / left wavefunctions
    // remember \int_0^\infty (e^{-\sqrt{2 E} x})^2 dx = 1/(2\sqrt{2 E})
    let norm_r_sqr = psi_r.iter()
        .map(|psi| psi.abs().powi(2))
        .sum::<Float>()
        * dx
        + 1./(2.*(-2.*E).sqrt());
    let norm_l_sqr = psi_l.iter()
        .map(|psi| psi.abs().powi(2))
        .sum::<Float>()
        * dx
        + 1./(2.*(-2.*E).sqrt());

    // We must now find the constants with which we should multiply
    // psi_l and psi_r with such that when we concatinate them,
    // we get a normalized and continuous wavefunction.
    //
    // psi = a * psi_l + b * psi_r
    // a / b = psi_r[-1] / psi_l[-1] ---- matching
    // a^2 * norm_l_sqr + b^2 * norm_r_sqr = 1 ---- normalized
    // => b = a * psi_l[-1] / psi_r[-1]
    // a = sqrt(1 / (norm_l_sqr + norm_r_sqr * psi_l[-1]^2 / psi_r[-1]^2 ))
    let a = -1. / (norm_l_sqr + norm_r_sqr * psi_l[psi_l.len()-1].powi(2) / psi_r[psi_r.len()-1].powi(2)).sqrt();
    let b = a * psi_l[psi_l.len()-1] / psi_r[psi_r.len()-1];


    // Calcualte the log derivative (for deciding if this is even a valid wf)
    let psi_l_prim_at_b = a * (psi_l[psi_l.len()-1] - psi_l[psi_l.len()-2]) / dx;
    let psi_r_prim_at_b = b * (psi_r[psi_r.len()-2] - psi_r[psi_r.len()-1]) / dx;
    let log_diff = (psi_l_prim_at_b - psi_r_prim_at_b) / psi_l[psi_l.len()-1];

    // Create the exponentially decreasing tails of the wavefunction
    // outside the potential
    let nbounds = (((support.0 - xbounds.0) / dx).floor() as usize,
        ((xbounds.1 - support.1) / dx).ceil() as usize);

    let psi_l_tail = (2..nbounds.0)
        .map(|n| dx * n as Float)
        .map(|x| (-(-2.*E).sqrt()*x).exp())
        .rev();
    let psi_r_tail = (2..nbounds.1)
        .map(|n| dx * n as Float)
        .map(|x| (-(-2.*E).sqrt()*x).exp());

    // Create the final wavefunction
    let psi_wf = psi_l_tail.map(|psi| psi*a)
        .chain(psi_l.into_iter().map(|psi| psi*a))
        .chain(psi_r.into_iter().rev().skip(1).map(|psi| psi*b))
        .chain(psi_r_tail.map(|psi| psi*b))
        .collect();

    let psi = State{
        energy: E,
        wf: psi_wf,
        x_left: xbounds.0,
        x_right: xbounds.1,
    };

    //println!("i: {}, length of pot: {}, energy_diff", i, potential.len());

    //println!("{:?}", psi.wf.iter().map(|psi| psi.abs().powi(2)).sum::<Float>()*dx);

    (psi, log_diff)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn numerov() {
        let potential = (0..100)
            .map(|i| 400. * ((i as Float / 100.).powi(2) - (i as Float / 100.)))
            .collect();
        let psi = super::numerov(-200., 0.99, 1., 1, 98, 0.001, &potential);
        // It should look smooth...
        println!("{:?}", psi);

    }
    #[test]
    fn bidir() {
        //let potential = (0..100)
            //.map(|i| 400. * ((i as Float / 100.).powi(2) - (i as Float / 100.)))
            //.collect();
        let potential = vec![-100.; 100];
        let psi = bidirectional_shooting(-50., (-1.0, 1.0), (-0.5,0.5), &potential);
        // It should look smooth...
        println!("{:?}", psi);
        //assert!(false);
    }

    #[test]
    fn find_bound_states_test() {
        let potential = vec![-100.; 100];
        let psi = find_bound_states_test((-1.0, 1.0), (-0.5,0.5), &potential);
        // It should look smooth...
        println!("{:?}", psi);
        assert!(false);
    }


}
