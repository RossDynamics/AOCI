function points = aeomega2r(a,e,omega,c)
%AEOMEGA2R Converts initial conditions whose a, e, and omega values are
%known and which must lie along the m1-l4 line to a position/momentum
%rotating frame state vector representation. The primary is assumed to be
%m_1. Vectors of a and e may be provided, but omega is constant for all
%initial conditions. A context object must be provided.

%we have to distinguish between the 3 body and 2 body mu here
mu3bdy = cg(c,'p.mu');
mu2bdy = 1 - mu3bdy;

%We calculate nu from omega, based on the fact that in the frame centered
%on the primary initial conditions are constrained to lie along 
%Y = sqrt(3) * X.
nufinder = @(nu) (sqrt(3)*(-cos(nu)*cos(omega) + sin(nu)*sin(omega)) ...
                + cos(nu)*sin(omega) + sin(nu)*cos(omega));
     
%Because of the polarity of the Poincare section, we assume v_r > 0.
%We use an initial guess of 4 to try to ensure the initial condition is in
%quadrant I w.r.t. the m_1-centered rotating frame coordinate system.
nu = fzero(nufinder,4);

%We get the location of the primary
xprime = - mu3bdy;

xlocs = [xprime*ones(1,size(a,2))
         zeros(1,size(a,2))];
     
%We have to compute h:
h = sqrt(mu2bdy .* (1-e.^2).*a);

%We now compute the state vectors in the perifocal frame:
rconst = h.^2./(mu2bdy .* (1 + e .* cos(nu)));
vconst = mu2bdy ./ h;
rvec = rconst .* [cos(nu)*ones(1,size(a,2)) 
                  sin(nu)*ones(1,size(a,2))];
vvec = vconst .* [-sin(nu)  *ones(1,size(a,2))
                 (e + cos(nu))*ones(1,size(a,2))];

%Now, we rotate the state vectors by omega, which corresponds to the
%transpose of the omega rotation matrix, in order to transform into the
%standard primary-centered inertial frame
R = [cos(omega) -sin(omega)
     sin(omega) cos(omega)];

rvecpri = R * rvec;
vvecpri = R * vvec;

%Now, we transform into the barycentric frame. Recall that the inertial 
%velocities coincide with the rotating frame momenta at t=0:
points = [rvecpri + xlocs
          vvecpri];

end

