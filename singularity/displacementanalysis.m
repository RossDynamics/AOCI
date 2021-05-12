%Attempts to analyze what happens when we make small displacements from the
%singularity.

%We create a circle of initial conditions that will be used to create the
%positions and velocities.
%Sun-Jupiter
%n = 1;
%pdelta = 1e-3;
%vdelta = 1.551669;
%displacements = onesphere(n);

n = 3;
%Let's construct a trajectory from the moon's surface. In nondimensional
%units, the distance from the center of the moon to its surface is the
%following pdelta:
pdelta = 0.00451974755594;
%We calculate the vdelta to use as the escape velocity:
vdelta = sqrt(2*cg(c,'p.mu')/pdelta);

%out of curiosity, we get the initial velocity, which will be the deltaV
%required for this launch, in km / s:
vdeltakms = vdelta * 1.025
displacements = onesphere(n);

conditions = [1 - cg(c,'p.mu') + pdelta * displacements(:,1)'
              pdelta * displacements(:,2)'
              vdelta * displacements'
              zeros(1,size(displacements,1))]
         
c = startCaching(c);

endtime = 4*2*pi;

hold on

for i = 1:size(conditions,2)
    integplot(linspace(0,endtime,100000),conditions(:,i),c)
end
axis equal

c = stopCaching(c);