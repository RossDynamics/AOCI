function c = LCR_PCR3BP_Context
%LCR_ER3BP_CONTEXT A constructor/profile for the PCR3BP as expressed using
%the Levi-Civita regularization.

c = Context(4);

syms u1(t) u2(t) udot1(t) udot2(t) U1(t) U2(t) t E mu
c = variableNames([u1 u2],[udot1 udot2],[U1 U2],t,c);

%We don't separate into kinetic and potential energy here... instead, we
%start from the Hamiltonian and stay entirely within momentum coordinates.
H = LCR_PCR3BP_H(u1,u2,U1,U2,E,mu);

c = cs(c,'d.H',simplify(H));

c = solveDynamicsFromH(c);

end