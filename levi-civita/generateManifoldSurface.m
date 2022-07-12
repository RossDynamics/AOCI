%Integrates the stable manifold to the m_2 singularity backwards to the
%surface of section used to generate the arches of chaos plots.

eHandle = cg(c,'d.Ehandle');

%The Sun-Jupiter distance in km.
L = 7.784e8;
%The length of an astronomical unit expressed in km.
AUinkm = 149597870.691;
%The Sun-Jupiter distance in AU.
%LinAU = L/AUinkm;
LinAU = 5.202850321292242;
%The radius of Jupiter in km.
Jrkm = 71492;
%The radius of Jupiter in denormalized units.
Jr = Jrkm / L;

T = 11.829163;
Tscale = T/(2*pi);

%This is used to put an upper bound on how long the integrator is allowed
%to integrate before it is interrupted. It should ideally be rather large. 
explorerEndTimeInLC = 0.2;

mu = cg(c,'p.mu');

%We create a circle of initial conditions in position space
numPos = 50;
posDisplacement = 1e-7;

thetas = linspace(0,2*pi,numPos);
dispCircle = [cos(thetas); sin(thetas)];

energy = -1.3;

c = cs(c,'p.E',energy);

c = startCaching(c);

%We call getFromCache for the sole purpose of caching the equations of
%motion.
[~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

surfpts = [];

times = linspace(0,explorerEndTimeInLC,50);
 
for i = 1:size(dispCircle,2)
    
    vec = dispCircle(:,i);
    
    positions = posDisplacement*vec + [1 - mu; 0];
    
    %We refine a velocity at the desired energy
    velAtEnergy = fzero(@(vel)(eHandle(0,...
        [positions;...
        vel*vec + [0 -1; 1 0]*positions;0])...
        -energy), 10);
    
    %         disp('Refined velocity:')
    %         disp(velAtEnergy)
    velocities = velAtEnergy*vec;
    
    manifoldICstandard = [positions;...
        velocities + [0 -1; 1 0]*positions];
    %We convert into Levi-Civita form
    [~,manifoldICLC] = standard2lc(0,manifoldICstandard,mu);
    
    %We now integrate the trajectory in the Levi-Civita frame
    try
        solbd = integ(times,[manifoldICLC;0],c);
    catch exception
        disp('An exception occurred while integrating. Ignoring.')
        disp('Exception details:')
        exception
        continue;
    end
    [ttrajbd,ytrajbd] = lc2standard(times,deval(solbd,times),mu);
    
    surfpts(i,:,:) = ytrajbd;

end

c = stopCaching(c);

xmat = reshape(surfpts(:,1,:),[],numel(times));
ymat = reshape(surfpts(:,2,:),[],numel(times));
zmat = reshape(surfpts(:,3,:),[],numel(times));

%surf(xmat,ymat,zmat,'FaceLighting','gouraud','FaceColor','interp')
