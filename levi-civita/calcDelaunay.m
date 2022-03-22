function [G,L,l,gbar] = calcDelaunay(points, primary, c)
%CALCDELAUNAYANGULAR Calculates Delaunay variables
%for the points provided with respect to the primary provided.

mu = cg(c,'p.mu');

%We get various orbital elements and representations of the
%position w.r.t. m_1:
[a,e,omega] = getae(points,primary,1,c);
polar = cart2polar(points,(1-primary) - mu);

r = polar(1,:);
theta = polar(2,:);
rdot = polar(3,:);
thetadot = polar(4,:);

G = r.^2.*(1+thetadot);
L = 1./sqrt(-(G.^2./r.^2)+2./r - rdot.^2);
l = acos((1-r./a)./e)-r.*rdot./L;
gbar = theta - acos((G.^2./r - 1)./e);

end

