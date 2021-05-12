function c = Sun_Jupiter_ER3BP_Context()
%SUN_JUPITER_ER3BP_CONTEXT Creates an ER3BP context for the Sun-Jupiter
%system.

c = ER3BP_Context;

c = cs(c,'p.mu',9.537e-4);
%We'll turn eccentricity on at some point.
c = cs(c,'p.e',0);

syms f(t) e
fEqn = diff(f,t) == (1 + e*cos(f))^2/((1-e^2)^(3/2));

c = addIntegratedParameter(f,fEqn,c);

end


