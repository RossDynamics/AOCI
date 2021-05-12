function c = Earth_Moon_ER3BP_e0_Context()
%EARTH_MOON_ER3BP_Context Creates an ER3BP context for the Earth-Moon
%system. This context has 0 eccentricity.

c = ER3BP_Context;

c = cs(c,'p.mu',0.0121505816234336);
c = cs(c,'p.e',0);

syms f(t) e
fEqn = diff(f,t) == (1 + e*cos(f))^2/((1-e^2)^(3/2));

c = addIntegratedParameter(f,fEqn,c);

end