function [r,theta,rdot,rthetadot]=xy2pol(XX);
% function [r,theta,rdot,rthetadot]=xy2pol(XX);
%
% Transforms from cartesian coordinates (x,y,xv,yv) to polar coordinates
% r,rdot,0<theta<2*pi,rthetadot(perpendicular velocity)
%
% Shane Ross (revised 6.24.97)
 
% Use only planar coords if XX, state vector matrix, is in 3D space
[m,N]=size(XX);
if     N==6, XX=[XX(:,1)            XX(:,2)   XX(:,4) XX(:,5)];
elseif N==3, XX=[XX(:,1) zeros(size(XX(:,1))) XX(:,2) XX(:,3)]; end

x=XX(:,1);y=XX(:,2);xv=XX(:,3);yv=XX(:,4);
v=[xv yv];

[theta,r]=cart2pol(x,y);
for k=1:length(theta),if theta(k)<0,theta(k)=theta(k)+2*pi;end;end;

rhatx=x./r;
rhaty=y./r;
rhat=[rhatx(:) rhaty(:)];

thetahatx=-y./r;
thetahaty= x./r;
thetahat=[thetahatx(:) thetahaty(:)];

rdot=v(:,1).*rhat(:,1) + v(:,2).*rhat(:,2);
rthetadot=v(:,1).*thetahat(:,1) + v(:,2).*thetahat(:,2);
