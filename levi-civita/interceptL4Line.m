function [position,isterminal,directions] = interceptL4Line(~,y)
%withinRadius An ode event function that detects crossings of the m1-L4 
%line, but in Levi-Civita coordinates.

%The line in (X,Y) coords is Y = sqrt(3) * X + sqrt(3), which can be
%converted easily into Levi-Civita coordinates to yield the expression
%below
position = sqrt(3)*y(1)^2 - sqrt(3)*y(2)^2 - 2*y(1)*y(2) + sqrt(3);

%We only go for interceptions that are
%on the side of the line closest to L4.
isterminal = y(1) * y(2) > 0;
directions = -1;
end

