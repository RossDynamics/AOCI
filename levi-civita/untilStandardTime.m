function [position,isterminal,directions] = untilStandardTime(t,y,endTime)
%withinRadius An ode event function that detects when the standard time
%endTime has been reached and then stops integration.
%Need to use SUN_JUPITER_ER3BP_CONTEXT_INTEGTIME or something similar that
%integrates the standard time alongside the equations of motion.

position = y(5) - endTime;

isterminal = 1;
directions = 0;
end

