function i = r2i(t,r)
%R2I Similar to Dr. Ross's rot2iner, but simpler and accompanied by i2r.
%Requires a 4D phase space expressed in column form, as opposed to a 6D
%phase space expressed in row form. t is the time; r are the coordinates.

rpos = r(1:2,:);
rvel = r(3:4,:);

A = [cos(t) -sin(t)
     sin(t)  cos(t)];

J = [0 1 
    -1 0];

ipos = A*rpos;
ivel = -A*J*rpos + A*rvel;

i = [ipos; ivel];

end

