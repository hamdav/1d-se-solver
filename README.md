Controls: 

- Draw while holding down the mouse button to draw a potential
- Q calculates the bound states for the potential 
	(and also fixes the potential so that its a valid function)
- C clears everything (except energy scaling)
- J and K increases / decreases the energy scaling respectively,
	effectively increasing / decreasing the depth of the potential.
	They also automatically trigger a recalculation of the bound states.
- Up / Down arrows show the above / below wavefunction

Problem:
The algorithm is unstable. When you draw something resembeling a double well and 
increase the energy scaling, eventually the algorithm fails. This is because
near the places where the nodes change, the log-diff diverges do infinity very
very very close to this place. Like $10^{-30}$ from it. I've tried moar
precision but that strategy barely even delayed the problems. I think what is
needed is a new solution. One that isn't integration. I think that finding a
better integration scheme won't help because the integration is _supposed_ to
diverge for energies that are not exactly the right energies. So better
integration schemes will probably just make sure even more that the energy needs
to be EXACTLY correct to get a normalizable wavefunction.

Solution?
Maybe completely rethink how I do the computation. For example, do exact
diagonalization and have the potential be an infinite potential well, with the
infinite walls being hidden out of sight. How bad could it be? If I cut it off
at the point where the screen contains 1000 wavelengths, I'd need to diagonalize
a 1000 x 1000 matrix. How hard is that? Probably doable.
