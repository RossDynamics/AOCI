function cart = m1polar2cart(polar, mu)
%M1POLAR2CART Converts initial conditions in polar coordinates (r, theta
%rdot, thetadot) centered at the m1 singularity to regular, barycentric,
%rotating frame cartesian momentum coordinates (x, y, px, py).

r = polar(1,:);
theta = polar(2,:);
rdot = polar(3,:);
thetadot = polar(4,:);

x = r.*cos(theta) - mu;
y = r.*sin(theta);
vx = rdot.*cos(theta) - r .* thetadot .* sin(theta);
vy = rdot.*sin(theta) + r .* thetadot .* cos(theta);

px = vx - y;
py = vy + x;

cart = [x y px py].';

end

