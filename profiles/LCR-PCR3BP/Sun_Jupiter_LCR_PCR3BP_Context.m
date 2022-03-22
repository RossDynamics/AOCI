function [c,c1] = Sun_Jupiter_LCR_PCR3BP_Context()
%SUN_JUPITER_ER3BP_CONTEXT Creates a PCR3BP context for the Sun-Jupiter
%system, expressed using the Levi-Civita regularization.

c = LCR_PCR3BP_Context;

c = cs(c,'p.mu',9.537e-4);

%We make a Sun-Jupiter context in the standard problem in order to retrieve
%the energy handle (in momentum coordinates!) for computing E.
c1 = Sun_Jupiter_ER3BP_Context;
c1 = useMomentum(c1);

%I *think* you need to set E for each trajectory.

c = useMomentum(c);
c = cs(c,'d.Ehandle',getEnergyHandle(c1));

end


