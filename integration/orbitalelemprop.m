function sol = orbitalelemprop(~,tspan,y0,~,numPoints,endRadius,mu,...
    primary,hitNum,feqn)
%ORBITALELEMPROP Propagates y0 in position/velocity coords by converting to
%orbital elements. Will stop integrating either when tspan(end) - tspan(1)
%is reached or the radius with respect to body crosses endRadius, whichever
%comes first. The solution struct will be evaluated over numPoints points.
%Be sure to convert to a function handle that integ can accept.
%Because this function ignores standard event handling, it is important to
%"tie" it to an event number using hitNum. It's also important to provide
%the equation of motion for integrating f.

%We convert the end time into a duration.
t1 = tspan(end) - tspan(1);

%We get the position and 2-body mu of the primary, which should expressed
%as either a 1 or a 2.
primex = (primary - 1) - mu;
mu2bdy = abs((2 - primary) - mu);

%Be sure y0 is already in position/velocity form! We create an
%(approximated as) inertial coordinate system centered on the desired
%primary; we subtract off the primary's inertial position and velocity:
primoffset = r2i(tspan(1),[primex 0 0 0].');

y0i = r2i(tspan(1),y0)-primoffset;

rvec = [y0i(1); y0i(2); 0];
vvec = [y0i(3); y0i(4); 0];
kvec = [0 0 1];

r = norm(rvec);
v = norm(vvec);

%We get the energy and the H vector:
energy = 0.5*v^2-mu2bdy/r;
Hvec = cross(rvec, vvec);

%We get the semi-major axis:
a = -mu2bdy/(2*energy);

evec = (cross(vvec,Hvec)-(mu2bdy/r)*rvec)/mu2bdy;
%We get the eccentricity:
e = norm(evec);

%We get the inclination:
i = acos(dot(kvec,Hvec)/norm(Hvec));

%nvec = cross(kvec,Hvec)/norm(cross(kvec,Hvec));

%We get the true anomaly while doing a quadrant check:
prenu = dot(rvec,evec)/(r*e);
if dot(rvec,vvec)>=0
    nu = acos(prenu);
else
    nu = 2*pi-acos(prenu);
end

%We determine the eccentric anomaly and the time past the time
%of perigee passage:
E = 2 * atan(sqrt((1-e)/(1+e))*tan(nu/2));
t0 = (E - e*sin(E))*sqrt(a^3/mu2bdy);

hitEnd = false;

%If the apoapsis becomes greater than or equal to endRadius at some point
%or if the eccentricity is greater than or equal to 1, we
%calculate the time past the time of perigee passage at which the orbit
%intersects with the sphere defined by endRadius. If it is less than t1,
%we must use it as the end time instead.
if e >= 1 || (1+e)*a >= endRadius
    %We calculate the true anomaly at the endRadius. We are only interested
    %in the intersection along the path from perigee to apogee, and as a
    %consequence we don't have to do a quadrant conversion.
    newnuposs = acos((a*(1-e^2)-endRadius)/(endRadius*e));
    
    %We now attempt to get the time remaining until the intersection
    Enewposs = 2 * atan(sqrt((1-e)/(1+e))*tan(newnuposs/2));
    t1candidate = sqrt(a^3/mu2bdy) * (Enewposs-e*sin(Enewposs));
    
    if t1candidate < t1+t0
        t1 = t1candidate-t0;
        hitEnd = true;
        %newnuEnd = newnuposs;
        %EnewEnd = Enewposs;
    end
    
end

%Now, we need to integrate the f parameter of the system between the start
%and end times to ensure it remains accurate
    
[~,f] = ode113(feqn,linspace(tspan(1),tspan(1)+t1,numPoints),y0(5));

y1 = zeros(size(numPoints,2)-1,5);

times = linspace(0,t1,numPoints);

for i = 2:numPoints
    ti = times(i);
    
    %We now solve for the eccentric anomaly after epoch;
    %We have to add t1 to t0:
    syms Esym
    Enew = double(vpasolve(sqrt(mu2bdy/a^3)*(ti+t0) == Esym - e*sin(Esym)));
    
    %Enew = fsolve(@(E)sqrt(mu2bdy/a^3)*(ti+t0) - E + e*sin(E),1i,...
    %              optimset('Display','off'));
    
    %We now get the true anomaly (E is likely complex so we need to get
    %the real part of newnu; the complex part should be negligible anyhow):
    newnu = real(2 * atan(sqrt((1+e)/(1-e))*tan(Enew/2)));
    
    
    %and now the radius and p:
    newr = a*(1-e^2)/(1+e*cos(newnu));
    p = a*(1-e^2);
    
    %If we got NaN, then most likely both the numerator and the denominator
    %of the newr expression were 0. This problem can occur if, for example,
    %you are trying to propagate a radial trajectory. To solve it, we 
    %     if isnan(newr)
    %         
    %     end
    
    %In the perifocal frame, we can now compute the position and velocity:
    rperi = [newr*cos(newnu); newr*sin(newnu); 0];
    vperi = sqrt(mu2bdy/p)*[-sin(newnu); (e + cos(newnu)); 0];
    
    %We now rotate back into the standard inertial frame:
    pvec = evec / e;
    wvec = Hvec / norm(Hvec);
    R = [pvec cross(wvec,pvec) wvec];
    
    newrvec = R*rperi;
    newvvec = R*vperi;
    
    primoffsetnew = r2i(tspan(1)+ti,[primex 0 0 0].');
    
    %Now, we convert back to the rotating frame:
    y1 = [y1 [i2r(tspan(1)+ti,[newrvec(1:2); newvvec(1:2)]+primoffsetnew)
               f(i)]];
    
end
%y1 = [y1nof; f(end)];

%We manually build out a solution struct:
sol = struct;
sol.solver = 'AOCI/orbitalelemprop';
sol.x = tspan(1)+times;
sol.y = [y0 y1];

if hitEnd
    sol.ie = hitNum;
    sol.xe = sol.x(end);
    sol.ye = sol.y(:,end);
end

end

