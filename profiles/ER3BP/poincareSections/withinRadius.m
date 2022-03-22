function [position,isterminal,directions] = withinRadius(t,y,xloc,...
                                                 distances,isterminal,...
                                                 directions)
%withinRadius An ode event function that detects proximity to the point
%(xloc,0) within certain distance(s) and going in certain direction(s).
%Whether the event is terminal can be specified. 
%Can handle multiple radii at once. Be sure that the origin of position
%space is set to (0,0).

onesArray = ones(numel(distances),1);

position = sqrt((y(1)-xloc)^2+(y(2))^2)*onesArray - distances;

end

