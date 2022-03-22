function sepAngle = calcSepAngle(r,thetadot,mu)
%CALCSEPANGLE Takes initial condition(s) for a trajectory in 
%secondary-centered polar coordinates and attempts to calculate the
%two-body separation angle 2*theta_infinity between the asymptotes, 
%assuming the trajectory is hyperbolic with respect to the secondary. ic
%must have r = rp, where rp is the closest approach distance, which implies
%rdot = 0. %We don't actually need theta or rdot for this calculation, as 
%theta does not appear in the equations and the rdot.^2 term in the 
%velocity equation is zero because rdot is assumed zero.

%For the secondary, the two-body and three-body mu coincide, so it is not
%necessary to draw a distinction between the two.

vsquared = r.^2.*thetadot.^2;

energy = vsquared./2 - mu./r;
a = -mu./(2.*energy);

sepAngle = 2.*acos(a/(r-a));

end

