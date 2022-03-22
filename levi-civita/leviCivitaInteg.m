function [y,lcsol,c] = leviCivitaInteg(tspan,y0,c)
%LEVICIVITAINTEG Uses the Levi-Civita regularization to integrate a
%trajectory provided in standard form. Transformations are done "under the
%hood" so that one doesn't have to touch LC coordinates directly. lcsol
%does, however, contain the original solution struct in Levi-Civita form if
%it is required. c needs to be a *Levi-Civita* context object.

%We convert the initial condition to Levi-Civita coordinates:
[t0lc,y0lc] = standard2lc(tspan(1),y0,cg(c,'p.mu'));

%We get the energy of the trajectory so that integration in Levi-Civita
%coordinates functions properly:
Ehandle = cg(c,'d.Ehandle');
c = cs(c,'p.E',Ehandle(t0lc,[y0lc; t0lc]));

%We integrate the trajectory in Levi-Civita coordinates:
[sols, c] = integplot(tspan,y0lc,c)

[regularStart, regularIC =

end

