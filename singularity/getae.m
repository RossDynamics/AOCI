function [a,e,omega] = getae(points,primary,L,c)
%GETAE Gets the a and e values with respect to primary (can be either 1 or
%2) of a collection of phase space points. L is the factor used to convert
%nondimensional coordinates back into dimensional coordinates. Now also
%gets omega, the argument of periapsis!

%we have to distinguish between the 3 body and 2 body mu here
mu3bdy = cg(c,'p.mu');
if (primary == 1) 
    mu2bdy = 1 - mu3bdy;
elseif (primary == 2)
    mu2bdy = mu3bdy;
%We experimentally allow calculation w.r.t. the barycenter
elseif (primary == 0)
    mu2bdy = 1;
end

%We get the location of the primary
xprime = (primary - 1) - mu3bdy;

if ~primary
    xprime = 0;
end

xlocs = [xprime*ones(1,size(points,2))
         zeros(1,size(points,2))];

%We get the r vector of each point relative to the primary and its
%magnitude
rvecs = [points(1:2,:) - xlocs
         zeros(1,size(points,2))];

r = vecnorm(rvecs);

%We get the velocity vectors and their magnitudes, but we have to calculate
%them in the inertial frame. We assume, without loss of generality, t=0,
%which interestingly enough means that the inertial velocities coincide
%with the rotating frame momenta
vvecs = [points(3,:)
         points(4,:)
         zeros(1,size(points,2))];

v = vecnorm(vvecs);

%We get the two-body energies:
energies = v.^2 ./ 2 - mu2bdy ./ r

%We get the semimajor axes in dimensional units:
%a = -mu2bdy./(2.*energies)./L;

%We get the semimajor axes in nondimensional units:
a = -mu2bdy./(2.*energies);

%We get the H vectors:
Hvecs = cross(rvecs,vvecs);

%We get the ecccentricity vectors and the eccentricities:
evecs = (cross(vvecs,Hvecs)-mu2bdy.*rvecs./r)./mu2bdy;

e = vecnorm(evecs);

omega = atan2(evecs(2,:),evecs(1,:));

%This little bit of logic accounts for whether the orbit is clockwise or
%counterclockwise
HzIsNegative = (Hvecs(3,:) < 0);
omega = 2*pi*HzIsNegative + (-1).^HzIsNegative .* omega;

end

