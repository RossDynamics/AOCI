function c = LCR_PCR3BP_Context_6D
%LCR_ER3BP_CONTEXT A constructor/profile for the PCR3BP as expressed using
%the Levi-Civita regularization that includes t and E as phase space variables.

c = Context(6);

syms u1(tau) u2(tau) t(tau) udot1(tau) udot2(tau) tdot(tau) U1(tau) ...
     U2(tau) E(tau) tau mu
c = variableNames([u1 u2 t],[udot1 udot2 tdot],[U1 U2 E],tau,c);

%We don't separate into kinetic and potential energy here... instead, we
%start from the Hamiltonian and stay entirely within momentum coordinates.
H = LCR_PCR3BP_H(u1,u2,U1,U2,E,mu);

c = cs(c,'d.H',simplify(H));

c = solveDynamicsFromH(c);

end