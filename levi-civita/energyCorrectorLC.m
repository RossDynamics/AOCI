function correctedIC = energyCorrector(r,theta,c)
%ENERGYCORRECTOR Builds, in Cartesian position/velocity coordinates,
%an IC from the provided r and theta coordinates 
%whose radial velocity corresponds to an energy of zero in
%the Levi-Civita regularization. Be sure to set the desired energy (when
%transforming out of LC coordinates) first. 
       
energyHandle = getEnergyHandle(c);

%We guess an initial radial velocity
radialGuess = 0.1;

function en = radialEnergy(radialVel)
    
    ineroffset = [r*cos(theta)
                  r*sin(theta)
                  radialVel*cos(theta)
                  radialVel*sin(theta)];

    en = energyHandle(0,ineroffset);
end

radialVel = fzero(@radialEnergy,radialGuess);

correctedIC = [r*cos(theta)
               r*sin(theta)
               radialVel*cos(theta)
               radialVel*sin(theta)];
end

