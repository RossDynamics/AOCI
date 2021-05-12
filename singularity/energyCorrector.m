function correctedIC = energyCorrector(r,theta,energy,mass,c)
%ENERGYCORRECTOR Builds, in Cartesian position/velocity coordinates,
%an IC from the provided r and theta coordinates about the specified mass
%(either 1 or 2) whose radial velocity corresponds to the desired energy.

mu = cg(c,'p.mu');

%A fun trick to avoid an if statement
bodypos = [(mass-1)-mu
           0];

position = [bodypos(1) + r*cos(theta)
            bodypos(2) + r*sin(theta)];

energyHandle = getEnergyHandle(c);

%We guess an initial radial velocity
radialGuess = 1;

function en = radialEnergy(radialVel)
    attempt = [position
                    radialVel * cos(theta)
                    radialVel * sin(theta)
                    0];
    en = energyHandle(0,attempt) - energy;
end

radialVel = fzero(@radialEnergy,radialGuess);

correctedIC = [position
               radialVel * cos(theta)
               radialVel * sin(theta)
               0];

end

