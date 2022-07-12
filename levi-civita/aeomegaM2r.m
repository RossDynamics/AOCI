function points = aeomegaM2r(a,e,omega,M,epochAngle,c)
%AEOMEGA2R Converts initial conditions whose a, e, omega, and M values are
%known to a position/momentum rotating frame state vector representation. 
%The primary is assumed to be m_1. Vectors of a and e may be provided, but
%omega and M are constant for all initial conditions. An epoch angle should
%also be provided to harmonize with the Rosengren group's work.
%A context object must be provided.

%we have to distinguish between the 3 body and 2 body mu here
mu3bdy = cg(c,'p.mu');
mu2bdy = 1 - mu3bdy;

%We calculate nu from M and e
nufinder = @(nu,e1)(atan2(sqrt(1-e1^2)*sin(nu)/(1+e1*cos(nu)),...
     (e1+cos(nu))/(1+e1*cos(nu)))-e1*sqrt(1-e1^2)*sin(nu)/(1+e1*cos(nu))-M); 
     
nu = zeros(size(e));

%We use an initial guess of 4 to try to ensure the initial condition is in
%quadrant I w.r.t. the m_1-centered rotating frame coordinate system.
for i = 1:numel(e)
    if isnan(e(i))
        nu = NaN;
        continue;
    end
    if e(i) >= 1
        nu = NaN;
        continue;
    end
    try
        nu(i) = fzero(@(nu1)nufinder(nu1,e(i)),2.8);%,optimset('Display','iter'));
    catch ME
        if strcmp(ME.identifier,'MATLAB:fzero:InvalidFunctionSupplied')
            nu = NaN;
        else
            disp(ME)
            rethrow(ME)
        end
    end
end

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
%omega rotation matrix, in order to transform into the
%standard primary-centered inertial frame
R = [cos(omega) -sin(omega)
     sin(omega) cos(omega)];

rvecpri = R * rvec;
vvecpri = R * vvec;

t = epochAngle;

A = [cos(t) -sin(t)
     sin(t)  cos(t)];

%Now, we transform into the barycentric frame. Recall that the inertial 
%velocities coincide with the rotating frame momenta at t=0:
points = [A zeros(2,2)
          zeros(2,2) A]*[rvecpri + xlocs
          vvecpri];
%points = [rvecpri + xlocs
%          vvecpri];
      
end

