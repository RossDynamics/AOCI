function [position,isterminal,direction] = xsection(t,y)
%U4SECTION An ode event function that detects crossings at the x-axis with
%negative vy. Confusingly, this function has been adapted to work in
%Levi-Civita coordinates, so we actually detect crossings of the y-axis in
%order to detect crossings at the y-axis in standard coordinates (see
%Szebehely (1967)).

position = y(1);
isterminal = 0;
direction = 0;

end

