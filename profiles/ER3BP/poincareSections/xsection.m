function [position,isterminal,direction] = xsection(t,y)
%U4SECTION An ode event function that detects crossings at the x-axis with
%negative vy.

position = y(2);
isterminal = 0;
direction = -1;

end

