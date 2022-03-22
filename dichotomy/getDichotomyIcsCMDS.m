function c = getDichotomyIcsCMDS(numConds,lowere,uppere,a,c)
%GETDICHOTOMYICSCMDS A wrapper to getDichotomyIcsCMDS that also converts
%out of the perifocal reference frame and stores the resultant initial
%conditions in a context object c. The two-body mu is automatically
%calculated from the three-body mu of the context object (make sure you're
%using a dynamical system in which the three-body mu exists and is at
%p.mu!), but a nondimensional a must be provided. You may need to add
%additional initial conditions to the ends of each vector if your model
%requires extended parameters.

mu = cg(c,'p.mu');

%We have to perform the ics computation in a basis that's rotated 60
%degrees in both position space and velocity space.
c = cs(c,'ac.basis',[cos(pi/3) -sin(pi/3) 0          0
                     sin(pi/3)  cos(pi/3) 0          0
                     0          0         cos(pi/3) -sin(pi/3)
                     0          0         sin(pi/3)  cos(pi/3)]);

c = cs(c,'ac.origin',[-mu 0 0 0].')                                      
                                     
mu2body = 1 - mu;

ics = getDichotomyInitialConditions(numConds,lowere, uppere,...
                                             a, mu);

c = cs(c,'aoci.di.ics',ics);

c = coordreset(c);

end

