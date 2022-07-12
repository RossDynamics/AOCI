function [c,c1] = Sun_Jupiter_LCR_PCR3BP_Context_Integtime()
%SUN_JUPITER_ER3BP_CONTEXT_INTEGTIME Creates a PCR3BP context for the 
%Sun-Jupiter system, expressed using the Levi-Civita regularization. 
%This is an extension of the standard SUN_JUPITER_ER3BP_CONTEXT context
%object that integrates the standard time (usually called t, but here
%denoted via the variable time) alongside the LCR equations of
%motion.

[c,c1] = Sun_Jupiter_LCR_PCR3BP_Context;

syms time(t) u1(t) u2(t)
c = addIntegratedParameter(time,diff(time,t)==u1^2+u2^2,c);

end


