function [tlc,ylc] = standard2lc(tstandard,ystandard,mu)
%STANDARD2LC Converts t and y arrays to Levi-Civita coordinates from standard
%coordinates. Be sure that coordinates in the y array utilize the vertical
%format.

x = ystandard(1,:);
y = ystandard(2,:);
px =ystandard(3,:);
py =ystandard(4,:);

X = x - (1 - mu);
Y = y;
PX = px;
PY = py - (1 - mu);

usqr = sqrt(X.^2+Y.^2);

tol = 1e-8;

%When we invert the transformation between Levi-Civita and standard
%coordinates, we find that we technically have four choices for u1.
%Two of those choices will always be either zero or will possess an
%imaginary part, so we discard those; the remaining two are either positive
%or negative real numbers. In general we choose the
%positive result.
%This transformation can fail in certain situations.

u1 = sqrt((X+usqr)./2);
%We have to handle these cases using alternate equations or the conversion
%fails.
if u1 == 0
    u2 = sqrt(-X);
else
    u2 = Y./(2.*u1);
end

if u2 == 0
    U1 = 2.*usqr.*(PY.*u2./u1 + PX)./(u1+u2.^2./u1);
    U2 = (2.*usqr.*PY - U1.*u2)./u1;
else
    U2 = 2.*usqr.*(PY.*u1./u2 - PX)./(u2+u1.^2./u2);
    U1 = (2.*usqr.*PY - U2.*u1)./u2;
end

ylc = [u1
       u2
       U1
       U2];

tlc = zeros(1,numel(tstandard));
         
for i = 2:numel(tlc)
    tlc(i) = trapz(tstandard(1:i),1./usqr(1:i));         
    %tstandard = tlc .* usqr;
end

end

