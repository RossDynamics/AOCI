function r = i2r(t,i)
%I2R Converts inertial coordinates to rotating coordinates.
%Requires a 4D phase space expressed in column form, as opposed to a 6D
%phase space expressed in row form. t is the time; r are the coordinates.

ipos = i(1:2,:);
ivel = i(3:4,:);

Ainv = [cos(t) sin(t)
       -sin(t)  cos(t)];

J = [0 1 
    -1 0];

rpos = Ainv*ipos;
rvel = J*Ainv*ipos + Ainv*ivel;

r = [rpos; rvel];

end

