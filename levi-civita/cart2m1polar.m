function polar = cart2m1polar(cart, mu)
%CART2M1POLAR Converts initial conditions in regular, barycentric,
%rotating frame cartesian momentum coordinates (x, y, px, py) to polar
%coordinates (r, theta, rdot, thetadot) centered at the m1 singularity.
%THIS FUNCTION LOOKS FLAWED; MOMENTUM AND VELOCITY SEEM CONFUSED. USE
%CART2POLAR INSTEAD

X = cart(1,:) - mu;
Y = cart(2,:);
VX = cart(3,:);
VY = cart(4,:);

r = vecnorm([X; Y]);
theta = atan2(Y,X);
rdot = (X .* VX + Y .* VY) ./ r;
thetadot = (rdot .* cos(theta) - VX)./(r .* sin(theta));

polar = [r; theta; rdot; thetadot];

end

