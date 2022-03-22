function H = LCR_PCR3BP_H(u1,u2,U1,U2,E,mu)
%LCR_PCR3BP_H The Hamiltonian for the PCR3BP as expressed in the Levi-
%Civita regularization. Refer to Paez and Guzzo (2020) at
%https://doi.org/10.1016/j.ijnonlinmec.2020.103417
%for the form of the regularization used.

usqd = u1^2 + u2^2;

line1 = (1/8)*(U1+2*usqd*u2)^2 + (1/8)*(U2-2*usqd*u1)^2;

line2 = -(1/2)*usqd^3 - mu - usqd * (E + (1-mu)^2/2);

line3 = -(1 - mu)*usqd*(1/sqrt(1 + 2*(u1^2 - u2^2) + usqd^2) ...
        + u1^2 - u2^2);

H = line1 + line2 + line3;

% mu1 = 1 - mu;
% mu2 = mu;
% 
% D = 4*(Q1^2+Q2^2);
% 
% f = Q1^2-Q2^2-mu1;
% g = 2*Q1*Q2;
% 
% r1 = sqrt((f-mu2)^2 + g^2);
% r2 = sqrt((f-mu1)^2 + g^2);
% 
% Vhat = mu1 / r1 + mu2 / r2; 
% 
% dfgdQ2 = 4*Q1^2*Q2 + 4*Q2^3 + 4*Q2*mu1;
% dfgdQ1 = 4*Q2^2*Q1 + 4*Q1^3 - 4*Q1*mu1;
% 
% H = D*T+(1/2)*(P1^2 + P2^2 + P1*dfgdQ2 - P2*dfgdQ1) - D*Vhat;
end
