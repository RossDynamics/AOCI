function polar = cart2polar(cart, xloc)
%CART2POLAR Converts initial conditions in regular, barycentric,
%rotating frame cartesian momentum coordinates (x, y, px, py) to polar
%coordinates (r, theta, rdot, thetadot) centered at the position space 
%point (xloc, 0).

x = cart(1,:);
y = cart(2,:);
vx = cart(3,:) + cart(2,:);
vy = cart(4,:) - cart(1,:);

X = x - xloc;
Y = y;
VX = vx;
VY = vy;

r = vecnorm([X; Y]);
theta = atan2(Y,X);
rdot = (X .* VX + Y .* VY) ./ r;
thetadot = -VX.*(Y./r.^2) + VY.*(X./r.^2);

polar = [r; theta; rdot; thetadot];

end

