function [position,isterminal,directions] = toAOCSection(t,y,mu,omega,M,epochAngle,Tscale,endTime)
%withinRadius An ode event function that intercepts the surface of section
%used to generate the initial conditions for the arches of chaos plots.
%endTime is a standard-time endTime after which to stop integrating.

isterminal = [1 1 1 1 1];
directions = [0 0 0 0 0];

%If a trajectory's projection into position space escapes past the
%sanityEscapeDistance from the barycenter, we stop integrating.
sanityEscapeDistance = 10;
%Similarly, if a trajectory gets within sanityEncounterDistance of the m1
%singularity, we stop integrating.
sanityEncounterDistance = 1e-3;
%If the trajectory goes under a certain threshold, it might try to wrap
%around resulting in a false hit, so we catch it before it does.
pos1Threshold = -5;

[~,yst] = lc2standard(t,y,mu);

mu2bdy = 1 - mu;
xlocs = [-mu;0];

%This is the location of the particle in a standard (non-LC) inertial
%reference frame
 
%The angle by which to rotate in order to account for the epoch and 
%everything else is based on multiple factors:
rotateAngle = epochAngle;% + y(5)/Tscale;

Ainv = [cos(rotateAngle) sin(rotateAngle)
       -sin(rotateAngle) cos(rotateAngle)];

%Now, we transform into the inertial frame. Recall that the inertial 
%velocities coincide with the rotating frame momenta at t=0:
points = [Ainv zeros(2,2)
         zeros(2,2) Ainv]*yst(1:4) - [xlocs; 0; 0];

%points = yst(1:4);

%We also have to rotate according to omega, but *after* centering at m1:
Rinv = [cos(omega) sin(omega)
       -sin(omega) cos(omega)];

%We get the r vector of each point relative to the primary and its
%magnitude
rvecs = [Rinv*points(1:2,:)
         0];

r = vecnorm(rvecs);

%We get the velocity vectors and their magnitudes, but we have to calculate
%them in the inertial frame. Because we are in momentum coordinates, we can
%write the following simplification:
vvecs = [Rinv*points(3:4,:)
         0];
     
v = vecnorm(vvecs);

energy = v.^2 ./ 2 - mu2bdy ./ r;

%We get the semimajor axes in nondimensional units:
a = -mu2bdy./(2.*energy);

%We get the H vectors:
Hvecs = cross(rvecs,vvecs);

%We get the ecccentricity vectors and the eccentricities:
evecs = (cross(vvecs,Hvecs)-mu2bdy.*rvecs./r)./mu2bdy;

e = vecnorm(evecs);
if e>=1
    position = [inf ...
                norm(yst(1:2)) - sanityEscapeDistance ...
                norm(yst(1:2) - xlocs) - sanityEncounterDistance ...
                abs(y(5)) - abs(endTime) ...
                inf];
    %If the eccentricity becomes greater than 1, we stop integrating
%     disp('Eccentricity > 1 @')
%     disp(e)
%     disp('LC time')
%     disp(t)
    return;
end

%nu = acos((abs(a*(1-e)^2)-r)/(e*r));
nu = acos(evecs'*rvecs/(e*r));

eDotrIsNegative = (rvecs'*vvecs < 0);
nu = 2*pi*eDotrIsNegative + (-1).^eDotrIsNegative .* nu;

%  disp('True anomaly:')
%  disp(nu)

%Mcurr = sqrt(mu2bdy/a^3)*y(5)
%scatter(y(5),Mcurr,'.')

position1threshold = 2.5e-2;

%omegaTarget = -0.688;
omegaTarget = 0;

omegaCalc = atan2(evecs(2,:),evecs(1,:));

% disp('Omega:')
% disp(omegaCalc)

position1 = abs(atan2(sqrt(1-e^2)*sin(nu)/(1+e*cos(nu)),...
     (e+cos(nu))/(1+e*cos(nu)))-e*sqrt(1-e^2)*sin(nu)/(1+e*cos(nu))-M) +...
     abs(omegaCalc - omegaTarget) - position1threshold;

%The surface of section corresponds to M reaching the designated value:
position = [position1,...
     norm(yst(1:2)) - sanityEscapeDistance,...
     norm(yst(1:2) - xlocs) - sanityEncounterDistance,...
     abs(y(5)) - abs(endTime),...
     inf];%position1 - pos1Threshold];

% scatter(y(5),position(1),'.')
% disp('pos(1):')
% disp(position(1))
% disp('t:')
% disp(y(5))
% hold on;
% drawnow;

%disp(position(1))

end