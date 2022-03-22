function correctedIC = energyCorrector(r,theta,energy,mass,c)
%ENERGYCORRECTOR Builds, in Cartesian position/velocity coordinates,
%an IC from the provided r and theta coordinates about the specified mass
%(either 1 or 2) whose radial velocity corresponds to the desired energy.

mu = cg(c,'p.mu');

%A fun trick to avoid an if statement
bodypos = [(mass-1)-mu
           0];
       
bodyiner = r2i(0,[bodypos; 0; 0]);
       
energyHandle = getEnergyHandle(c);

%We guess an initial radial velocity
radialGuess = 14;

function en = radialEnergy(radialVel)
    
    ineroffset = [r*cos(theta)
                  r*sin(theta)
                  radialVel*cos(theta)
                  radialVel*sin(theta)];
              
   attempt = [i2r(0,bodyiner+ineroffset)
              0]; 
              
   % attempt = [position
   %                 radialVel * cos(theta)
   %                 radialVel * sin(theta)
   %                 0];
    en = energyHandle(0,attempt) - energy;
end

%radialVel = radialGuess;

radialVel = fzero(@radialEnergy,radialGuess);

ineroffset = [r*cos(theta)
               r*sin(theta)
               radialVel*cos(theta)
               radialVel*sin(theta)];
           
correctedIC = [i2r(0,bodyiner+ineroffset)
              0]; 

% correctedIC = [position
%                radialVel * cos(theta)
%                radialVel * sin(theta)
%                0];

end

