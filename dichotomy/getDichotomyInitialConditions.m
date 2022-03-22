function ics = getDichotomyInitialConditions(numConds,lowere, uppere,...
                                             a, mu)
%GETDICHOTOMYINITIALCONDITIONS Generates initial conditions along the
%perihelions of two-body orbits with the specified range of eccentricities
%and with the provided semi-major axis. Don't call this function directly
%as it doesn't convert out of the perifocal coordinate system. Instead use
%getDichotomyIcsCMDS, which leverages a context object's automatic
%coordinate conversion capabilities. mu should be the *two-body* mass
%parameter; it should be possible to just use the nondimensional mass of
%m_1 if nondimensional a is also provided.

eArray = linspace(lowere,uppere,numConds);

%We generate initial conditions at periapsis.
rArray = a .* (1 - eArray);

%An extremely minor optimization.
pArray = rArray .* (1 - eArray);

ics = [rArray
       zeros(1,numel(rArray))
       zeros(1,numel(rArray))
       sqrt(mu./pArray).*(1+eArray)];


end

