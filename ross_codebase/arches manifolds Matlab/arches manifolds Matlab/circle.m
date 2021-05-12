function circle(r,c,color,FILL,LineWidth)
% circle(r,c,color,FILL,LineWidth);
%
% Make circle of with specified 'color' of radius r 
% at center c=[x y] for a 2D plot
%           c=[x y z] for a 3D plot
%
% if FILL=1, the fill circle

if nargin<=4, LineWidth=1; end
if nargin<=3, FILL=0; end
if nargin<=2, color ='w:';if nargin==1,c=[0 0];end; end
 
theta = 0:.0025:2*pi; theta = theta(:);

m=length(c); for k = 1:length(theta), C(k,1:m)=c; end 

if FILL==0,
    if     m==2, circ = plot (C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),color);
    elseif m==3, circ = plot3(C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),C(:,3),color);
    LineWidth
    set(circ,'LineWidth',LineWidth)
    end
elseif FILL==1,
    fill(C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),color);
end
