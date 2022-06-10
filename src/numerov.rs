

#[derive(Debug)]
pub struct State{
    pub energy: f64,
    pub wf: Vec<f64>,
    pub x_left: f64,
    pub x_right: f64,
}

pub fn find_bound_states(xbounds: (f64, f64),
                     support: (f64, f64), 
                     potential: &Vec<f64>) -> Vec<State> 
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
    
    let E_min = *potential.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let energy_interval = (E_min, E_min/1000.);
    println!("E_min: {}, energy_interval: {:?}", E_min, energy_interval);
    let (psi_lo, log_diff_lo) = bidirectional_shooting(E_min, xbounds, support, potential);
    let (psi_hi, log_diff_hi) = bidirectional_shooting(E_min/1000., xbounds, support, potential);
    //println!("{:?}\nnodes: {}", psi_hi, count_nodes
    //
    let mut rv = Vec::new();

    for n_nodes in 0.. {
        let E = recursive_bisect(energy_interval, (log_diff_lo, log_diff_hi),
            (count_nodes(&psi_lo), count_nodes(&psi_hi)), n_nodes, xbounds, support, potential);
        if let Some(E) = E {
            let (psi, _) = bidirectional_shooting(E, xbounds, support, potential);
            rv.push(psi);
        } else {
            break;
        }
    }
    rv
} 

pub fn find_bound_states2(xbounds: (f64, f64),
    support: (f64, f64), 
    potential: &Vec<f64>) -> Vec<State> 
{
    let E_min = *potential.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let ediffn = (0..1000).map(|n| {
        let E = E_min * (1000. - n as f64) / 1000.;
        let (psi, log_diff) = bidirectional_shooting(E, xbounds, support, potential);
        (E, log_diff, count_nodes(&psi))
    })
    .collect::<Vec<(f64, f64, usize)>>();

    println!("{:?}", ediffn.iter().map(|(_,b,_)| *b).collect::<Vec<f64>>());

    let psis = ediffn.iter()
        .zip(ediffn.iter().skip(1))
        // Filter out only the crossings
        .filter(|((_, log_diff1, n_nodes1), (_, log_diff2, n_nodes2))|
                n_nodes1 == n_nodes2 && log_diff1 * log_diff2 <= 0.)
        // Find crossings with bisection method
        .map(|((E1, log_diff1, _), (E2, log_diff2, _))|
             simple_logdiff_bisect((*E1, *E2), (*log_diff1, *log_diff2),
                                   xbounds, support, potential))
        // Now we have an energy, make it a psi
        .filter_map(|opt_e| opt_e.map(|energy| bidirectional_shooting(energy, xbounds, support, potential).0))
        .collect();

    psis
}

pub fn find_bound_states3(xbounds: (f64, f64),
                     support: (f64, f64), 
                     potential: &Vec<f64>) -> Vec<State> 
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
    
    let E_min = *potential.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let energy_interval = (E_min, E_min/1000.);
    println!("E_min: {}, energy_interval: {:?}", E_min, energy_interval);
    let (psi_lo, log_diff_lo) = bidirectional_shooting(E_min, xbounds, support, potential);
    let (psi_hi, log_diff_hi) = bidirectional_shooting(E_min/1000., xbounds, support, potential);
    //println!("{:?}\nnodes: {}", psi_hi, count_nodes
    //
    let mut rv = Vec::new();

    let mut lower_end = E_min;
    for n_nodes in 0.. {
        let interval = find_interval_enclosing_upper_boundry(energy_interval,
            (count_nodes(&psi_lo), count_nodes(&psi_hi)), n_nodes, xbounds, support, potential);
        if let Some(interval) = interval {
            let (psi_lo, log_diff_lo) = bidirectional_shooting(lower_end, xbounds, support, potential);
            let (psi_hi, log_diff_hi) = bidirectional_shooting(interval.0, xbounds, support, potential);

            println!("");
            println!("E: {:?}", (lower_end, interval.0));
            println!("log_diffs: {:?}", (log_diff_lo, log_diff_hi));
            println!("# of nodes: {:?}", (count_nodes(&psi_lo), count_nodes(&psi_hi)));
            let E = simple_logdiff_bisect((lower_end, interval.0),
                (log_diff_lo, log_diff_hi),
                xbounds, support, potential);

            if let Some(E) = E {
                let (psi, _) = bidirectional_shooting(E, xbounds, support, potential);
                rv.push(psi);
            } else {
                println!("Failed logdiff bisection");
            }

            // Prepare for next loop
            lower_end = interval.1;
        } else {
            // Do one with the upper end being E_min / 1000.
            let (psi_lo, log_diff_lo) = bidirectional_shooting(lower_end, xbounds, support, potential);
            let (psi_hi, log_diff_hi) = bidirectional_shooting(energy_interval.1, xbounds, support, potential);

            println!("");
            println!("E: {:?}", (lower_end, energy_interval.1));
            println!("log_diffs: {:?}", (log_diff_lo, log_diff_hi));
            println!("# of nodes: {:?}", (count_nodes(&psi_lo), count_nodes(&psi_hi)));
            let E = simple_logdiff_bisect((lower_end, energy_interval.1),
                (log_diff_lo, log_diff_hi),
                xbounds, support, potential);

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

fn simple_logdiff_bisect(mut energy_interval: (f64, f64),
                         mut end_log_diffs: (f64, f64),
                         xbounds: (f64, f64),
                         support: (f64, f64),
                         potential: &Vec<f64>)
    -> Option<f64>
{
    /*
     * Returns Some(energy) if the different ends of the energy interval
     * has the same number of nodes and different signs for the log diff
     */

    // Check that interval is valid
    if end_log_diffs.0 * end_log_diffs.1 > 0. {
        return None
    }
    // Base case: interval is small = return middle of interval.
    let mid = (energy_interval.0 + energy_interval.1) / 2.;
    if ((energy_interval.1 - energy_interval.0) / energy_interval.0).abs() < 1e-14 {
        return Some(mid)
    }

    // Check mid
    let (_, log_diff) = bidirectional_shooting(mid, xbounds, support, potential);

    if log_diff * end_log_diffs.0 <= 0. {
        end_log_diffs.1 = log_diff;
        energy_interval.1 = mid;
    } else {
        end_log_diffs.0 = log_diff;
        energy_interval.0 = mid;
    }
    return simple_logdiff_bisect(energy_interval, end_log_diffs, xbounds, support, potential);
}

fn find_interval_enclosing_upper_boundry(mut energy_interval: (f64, f64),
                                   mut end_nodes: (usize, usize),
                                   n_nodes: usize,
                                   xbounds: (f64, f64),
                                   support: (f64, f64),
                                   potential: &Vec<f64>)
    -> Option<(f64, f64)>
{
    /*
     * Finds the upper endpoint of the energy interval where the wf
     * has n_nodes nodes
     */

    // Base case: interval is small = return lower end of interval.
    if ((energy_interval.1 - energy_interval.0) / energy_interval.0).abs() < 1e-14 {
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

    if mid_n_nodes <= n_nodes {
        end_nodes.0 = mid_n_nodes;
        energy_interval.0 = mid;
    } else {
        end_nodes.1 = mid_n_nodes;
        energy_interval.1 = mid;
    }

    return find_interval_enclosing_upper_boundry(energy_interval, end_nodes, n_nodes, xbounds, support, potential)
}

fn recursive_bisect(mut energy_interval: (f64, f64), mut end_log_diffs: (f64, f64),
                    mut end_nodes: (usize, usize),
                    n_nodes: usize, xbounds: (f64, f64), 
                    support: (f64, f64),
                    potential: &Vec<f64>)
    -> Option<f64>
{
    /*
     * recursively bisect to find the 0 of the bidir log diff stuff
     */
    // Base case: interval is small = return middle of interval.
    if ((energy_interval.1 - energy_interval.0) / energy_interval.0).abs() < 1e-5 {
        if end_nodes.0 == n_nodes && end_nodes.1 == n_nodes {
            return Some((energy_interval.0 + energy_interval.1) / 2.)
        } else {
            return None
        }
    }

    // This will keep track of which endpoint we should replace
    let mut replace_lo = false;

    let mid = (energy_interval.0 + energy_interval.1) / 2.;
    let (psi, log_diff) = bidirectional_shooting(mid, xbounds, support, potential);
    let mid_n_nodes = count_nodes(&psi);
    println!("");
    println!("E: {:?}", energy_interval);
    println!("log_diffs: {:?}", end_log_diffs);
    println!("# of nodes: {:?}", end_nodes);
    println!("mid: {}, {}, {}", mid, mid_n_nodes, log_diff);

    // Check if the bound has too many or to few nodes and cut it accordingly
    if mid_n_nodes < n_nodes {
        replace_lo = true;
    } else if mid_n_nodes > n_nodes {
        replace_lo = false;
    }

    // If mid and the lower end has the correct number of nodes,
    // the root is between them is they have opposite signs, otherwise
    // it's between mid and high end
    else if mid_n_nodes == n_nodes && end_nodes.0 == n_nodes {
        if log_diff * end_log_diffs.0 <= 0. {
            replace_lo = false;
        } else {
            replace_lo = true;
        }
    }
    // Same with mid and the upper end
    else if mid_n_nodes == n_nodes && end_nodes.1 == n_nodes {
        if log_diff * end_log_diffs.1 <= 0. {
            replace_lo = true;
        } else {
            replace_lo = false;
        }
    }


    // So, now we know mid_n_nodes == n_nodes but neither of the ends 
    // have the correct number of nodes.
    // We still don't know if the root is to the right or to the left though
    // because we're not sure if the function is monotonously increasing or 
    // decreasing.
    // So let's find out by just checking
    else {
        let mut trial_E = (energy_interval.1 + mid) / 2.;
        let mut func_is_inc = false;
        println!("Looping");
        for i in 0..25 {
        //loop {
            let (trial_psi, trial_log_diff) = bidirectional_shooting(trial_E, xbounds, support, potential);
            let trial_nodes = count_nodes(&trial_psi);
            println!("trial: {}, {}, {}", trial_E, trial_nodes, trial_log_diff);
            if trial_nodes == n_nodes {
                if trial_log_diff > log_diff {
                    func_is_inc = true;
                } else {
                    func_is_inc = false;
                }
                break;
            } else {
                trial_E = (trial_E + mid) / 2.;
            }
            if i == 24 {
                println!("{:?}", trial_psi);
                println!("{:?}", psi);
                panic!("öajsdöflkasdf");
            }
        }
        
        if func_is_inc {
            replace_lo = log_diff <= 0.;
        } else {
            replace_lo = log_diff >= 0.;
        }
    }


    if replace_lo {
        end_nodes.0 = mid_n_nodes;
        end_log_diffs.0 = log_diff;
        energy_interval.0 = mid;
        return recursive_bisect(energy_interval, end_log_diffs, end_nodes, n_nodes, xbounds, support, potential)
    } else {
        end_nodes.1 = mid_n_nodes;
        end_log_diffs.1 = log_diff;
        energy_interval.1 = mid;
        return recursive_bisect(energy_interval, end_log_diffs, end_nodes, n_nodes, xbounds, support, potential)
    }
}


fn min_index<I>(iter: I) -> usize
    where I: IntoIterator<Item = f64>
{
    let mut min = f64::INFINITY;
    let mut mindex = 0;

    for (i, v) in iter.into_iter().enumerate() {
        if v < min {
            min = v;
            mindex = i;
        }
    }
    mindex
}
                
fn numerov(E: f64, psi_0: f64, psi_1: f64, start_ind: usize, end_ind: usize, dx: f64,
           potential: &Vec<f64>) -> Vec<f64>
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
    let mut psis = Vec::with_capacity((end_ind as i32 - start_ind as i32).abs() as usize);
    psis.push(psi_0);
    psis.push(psi_1);

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
    


fn bidirectional_shooting(E: f64, 
                          xbounds: (f64, f64), 
                          support: (f64, f64), 
                          potential: &Vec<f64>)
    -> (State, f64)
{
    /*
     * runs numerov from left and from right to b
     */

    assert!(E < 0.);

    let dx = (support.1 - support.0) / (potential.len() - 1) as f64;

    // Find b as the point where the potential is closest to E
    // What if there is a discontinuity?
    // What if there are multiple degenerate wavefunctions at this energy?
    //let mut i = min_index(potential.iter().map(|v| (v-E).abs()));
    let mut i = potential.iter().take_while(|&&v| v > E).count();
    if i <= 1 { i = 2; }
    else if i >= potential.len() - 2 {i = potential.len()-3;}
    //i = potential.len() / 2;

    // if V[i] - E is too large, do something...
    let psi_l = numerov(E, (-(-2.*E).sqrt() * dx).exp(), 1., 1, i, dx, &potential);
    let psi_r = numerov(E, (-(-2.*E).sqrt() * dx).exp(), 1., potential.len()-2, i, dx, &potential);

    // Find the norms of the right / left wavefunctions
    // remember \int_0^\infty e^{-\sqrt{2 E} x} dx = 1/\sqrt{2 E}
    let norm_r_sqr = psi_r.iter()
        .map(|psi| psi.abs().powi(2))
        .sum::<f64>()
        * dx
        + 1./(2.*(-2.*E).sqrt());
    let norm_l_sqr = psi_l.iter()
        .map(|psi| psi.abs().powi(2))
        .sum::<f64>()
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


    // Calcualte the log derivative (for deciding if this is even a valid wf
    let psi_l_prim_at_b = a * (psi_l[psi_l.len()-1] - psi_l[psi_l.len()-2]) / dx;
    let psi_r_prim_at_b = b * (psi_r[psi_r.len()-2] - psi_r[psi_r.len()-1]) / dx;
    let log_diff = (psi_l_prim_at_b - psi_r_prim_at_b) / psi_l[psi_l.len()-1];

    let nbounds = (((support.0 - xbounds.0) / dx).floor() as usize,
        ((xbounds.1 - support.1) / dx).ceil() as usize);

    let psi_l_tail = (2..nbounds.0)
        .map(|n| dx * n as f64)
        .map(|x| (-(-2.*E).sqrt()*x).exp())
        .rev();
    let psi_r_tail = (2..nbounds.1)
        .map(|n| dx * n as f64)
        .map(|x| (-(-2.*E).sqrt()*x).exp())
        .peekable();

    //println!("psi_r: {}, psi_r_tail: {}", psi_r[0], *psi_r_tail.peek().unwrap());

    //println!("Psi_L: {:?}", psi_l);
    //println!("Psi_L: {:?}", psi_l);
    

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

    //println!("{:?}", psi.wf.iter().map(|psi| psi.abs().powi(2)).sum::<f64>()*dx);

    (psi, log_diff)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn numerov() {
        let potential = (0..100)
            .map(|i| 400. * ((i as f64 / 100.).powi(2) - (i as f64 / 100.)))
            .collect();
        let psi = super::numerov(-200., 0.99, 1., 1, 98, 0.001, &potential);
        // It should look smooth...
        println!("{:?}", psi);

    }
    #[test]
    fn bidir() {
        //let potential = (0..100)
            //.map(|i| 400. * ((i as f64 / 100.).powi(2) - (i as f64 / 100.)))
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
