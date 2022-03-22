function [tstandard,ystandard] = lc2standard(tlc,ylc,mu)
%LC2STANDARD Converts t and y arrays in Levi-Civita coordinates to standard
%coordinates. Be sure that coordinates in the y array utilize the vertical
%format.

u1 = ylc(1,:);
u2 = ylc(2,:);
U1 = ylc(3,:);
U2 = ylc(4,:);

usqr = u1.^2 + u2.^2;

X = u1.^2 - u2.^2;
Y = 2.*u1.*u2;
PX = (U1.*u1 - U2.*u2)./(2.*usqr);
PY = (U1.*u2 + U2.*u1)./(2.*usqr);

ystandard = [X + 1 - mu
             Y
             PX
             PY + 1 - mu];       
         
if numel(tlc) > 1        
    tstandard = cumtrapz(tlc,usqr); 
else
    tstandard = 0;
end
end

