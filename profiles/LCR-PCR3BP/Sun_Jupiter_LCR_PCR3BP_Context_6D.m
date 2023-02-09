function [c,c1] = Sun_Jupiter_LCR_PCR3BP_Context_6D()
%SUN_JUPITER_ER3BP_CONTEXT_6D Creates a PCR3BP context for the 
%Sun-Jupiter system, expressed using the Levi-Civita regularization. 
%This is an extension of the standard SUN_JUPITER_ER3BP_CONTEXT context
%object that integrates the standard time (usually called t, but here
%denoted via the variable time)*and* the standard energy E alongside the 
%LCR equations of motion. This context object is useful for analyzing phase
%space rigorously.

% [c,c1] = Sun_Jupiter_LCR_PCR3BP_Context;
% 
% syms time(t) u1(t) u2(t) E(t)
% c = addIntegratedParameter(time,diff(time,t)==u1^2+u2^2,c);
% c = addIntegratedParameter(E,diff(E,t)==0,c);

c = LCR_PCR3BP_Context_6D;

c = cs(c,'p.mu',9.537e-4);


end


