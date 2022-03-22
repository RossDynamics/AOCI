%A tester to ensure that the Levi Civita context is set up properly.

close all;

%c = Sun_Jupiter_LCR_PCR3BP_Context;
%c1 = Sun_Jupiter_ER3BP_Context;

%We provide the initial condition in Levi-Civita coordinates:
ic = [0.1 0.2 0.3 0.4].'
endTime = -4;
numPts = 1000;
%We convert the initial condition to standard coordinates:
[~,ystr] = lc2standard(0,ic,cg(c,'p.mu'))

% ystr = [1 - cg(c,'p.mu') + 0.1 0 0.1 0].'
% [~,ic] = standard2lc(0,ystr,cg(c,'p.mu'))

%We get the energy of the trajectory so that integration in Levi-Civita
%coordinates functions properly:
Ehandle = cg(c,'d.Ehandle');
c = cs(c,'p.E',Ehandle(0,[ystr; 0]));

%We integrate and plot the trajectory in Levi-Civita coordinates:
[p,c,lcsol] = integplot(linspace(0, endTime, numPts),ic,c);

%We convert the trajectory to standard coordinates and compare it to the
%regular way of integrating:
figure;
hold on;
[tstrtraj,ystrtraj] = lc2standard(lcsol.x,lcsol.y,cg(c,'p.mu'));
cplot(ystrtraj,c)
drawnow;

c1 = useMomentum(c1);
[~,~,standsol] = integplot(linspace(tstrtraj(1), tstrtraj(end), numPts),...
                           [ystr;0],c1);