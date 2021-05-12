function [q1,q2,p1,p2,a]=rot2delB(XX);
% function [q1,q2,p1,p2,a]=rot2delB(XX);
%
% Transform from CR3BP rotating coords to Delaunay (canonical) coords 
%
% input: XX=state vector matrix (rotating coords)
%
% output: q1=mean anomaly
%	  q2=argument of periapse in rotating frame
%	  p1=sqrt(a)
%	  p2=sqrt(a*(1-e^2))
%
%	where (a,e) are the semimajor axis and eccentricity of an ellipse
%	representing the motion of the particle when acted upon only by a unit 
%	mass at the center of mass (i.e., mu=mu2=0, mu1=1)
%
%	see p.363 of Szebehely
%
% Shane Ross (revised 9.05.01)
[r,theta,rdot,rthetadot]=xy2pol(XX); % polar coords in rotating frame
%for k=1:length(theta), if theta(k)>pi, theta(k)=theta(k)-2*pi;end;end

p2=r.*r + r.*rthetadot; % p2=r^2(1+thetadot)
%p2=r.*rthetadot; % p2=r^2*thetadot
a=1./( -((p2.*p2)./(r.*r)) + 2./r - rdot.*rdot ); 

p1=sqrt(a) ;
e=sqrt(1-(p2.^2)./a);

q1=acos((1-r./a)./e) - r.*rdot./p1;
q2=theta - acos((a.*(1-e.^2)./r - 1)./e);

% for theta=0 and XX=[x xdot ydot], 
%we use the following to the angle right
l=size(XX,2) ;
if l==3, m=2;
elseif l==4, m=3;
elseif l==6, m=4;
end
q1=-sign(XX(:,m)).*q1 ;
q2=-sign(XX(:,m)).*q2 ;

%q1=arccos((1-r./a)./e) - r.*rdot./p1;
%q2=theta - arccos((a.*(1-e.^2)./r - 1)./e);

z2pi = 0 ;
if z2pi == 1,
 ztol = 0 ;
 for k=1:length(q2), 
  if q2(k)<ztol, 
    q2(k)=q2(k)+2*pi;
  end
  if q1(k)<ztol, 
    q1(k)=q1(k)+2*pi;
  end
 end
end
%for k=1:length(q2)	% gives smooth angle as a function of time
%	if abs(q2(k))>pi, q2(k)=q2(k)-sign(q2(k))*2*pi; end
%	if q2(k)+pi<0,q2(k)=q2(k)+2*pi; end % make sure all within +-pi
%	if q1(k)+pi<0,q1(k)=q1(k)+2*pi; end % ditto
%end
