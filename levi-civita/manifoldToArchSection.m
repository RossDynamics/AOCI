%Integrates the stable manifold to the m_2 singularity backwards to the
%surface of section used to generate the arches of chaos plots.

aocomega = deg2rad(2.738332953159384e2);
M = deg2rad(10.684213694231330e1);
%M = deg2rad(4.684213694231330E+01);

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

%Our epoch is 30 September 2012
epochinyears = years(datetime(2012,09,30)-datetime(2000,1,1,12,0,0));
epochdenorm = epochinyears/Tscale;

%The J2000 mean longitude angle:
j2000Angle = deg2rad(34.40438);

%epochAngle = epochdenorm + j2000Angle;
epochAngle = 0.8159069576305198;

%This is used to put an upper bound on how long the integrator is allowed
%to integrate before it is interrupted. It should ideally be rather large. 
explorerEndTimeInLC = -500;
%This is the actual end time in nondimensional time units.
standardEndTime = -33 / Tscale

mu = cg(c,'p.mu');

%We create a circle of initial conditions in position space and a range of
%energies corresponding to each point along the circle from which the
%velocities will be generated
numPos = 200;
posDisplacement = 1e-7;

thetas = linspace(0,2*pi,numPos+1);
thetas = thetas(1:end-1);
dispCircle = [cos(thetas); sin(thetas)];

energyRange = linspace(-1.5194,0,800);
%energyRange = linspace(-1.1413,0,200);

c = cs(c,'s.i.odeopts',odeset('Events',@(t,y)toAOCsection(t,y,mu,aocomega,...
                      M,epochAngle,Tscale,standardEndTime),'RelTol',3e-14,'AbsTol',1e-15));
          
hitpts = zeros(numel(energyRange),size(dispCircle,2),4);
 
%We build out the initial conditions and integrate, choosing one energy at
%a time for speed
for j = 1:numel(energyRange)
    
    clf;

    energy = energyRange(j);

    c = cs(c,'p.E',energy);

    c = startCaching(c);

    %We call getFromCache for the sole purpose of caching the equations of
    %motion.
    [~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

    disp('Energy:')
    disp(energy)
    
    parfor i = 1:size(dispCircle,2)
        vec = dispCircle(:,i);
        
        positions = posDisplacement*vec + [1 - mu; 0];
        
        %We refine a velocity at the desired energy
        velAtEnergy = fzero(@(vel)(eHandle(0,...
              [positions;...
              vel*vec + [0 -1; 1 0]*positions;0])...
              -energy), -10);
          
%         disp('Refined velocity:')
%         disp(velAtEnergy)
        velocities = velAtEnergy*vec;
        
        
        manifoldICstandard = [positions;...
                              velocities + [0 -1; 1 0]*positions];
        %We convert into Levi-Civita form
        [~,manifoldICLC] = standard2lc(0,manifoldICstandard,mu);
        
        %We now integrate the trajectory in the Levi-Civita frame
        try
            solbd = integ(linspace(0,explorerEndTimeInLC,1000),[manifoldICLC;0],c);
        catch exception
            disp('An exception occurred while integrating. Ignoring.')
            disp('Exception details:')
            exception
            continue;
        end
        [ttrajbd,ytrajbd] = lc2standard(solbd.x,solbd.y,mu);
        
        %We permit NaN trajectories to pass silently
%         if ~numel(solbd.ie) && ~any(isnan(solbd.y),'all')
%             throw(MException('AOCI:ExplorerEndTimeInLCTooSmall',...
% 'explorerEndTimeInLC needs to be larger so that integration can finish.'));
%         end
        
        %endtimestandard(j,i) = ttrajbd(end);
        
        %Now, we collect the hits to the surface of section. 
        %We throw out trajectories that fail due to a variety of conditions
        if ~isempty(solbd.xe) && ~any(ismember([2 3 4 5],solbd.ie))
            hitpts(j,i,:) = ytrajbd(:,end);
            disp('Hit')
            disp(i)
            disp(solbd.y(5,end))
            
%             cplot(ytrajbd,c);
%             cplot(ytrajbd(:,end),c,'o');
%             hold on;
%             drawnow;

        end
        
        disp(solbd.ie)
        
         
    end
    
    flathitpts = reshape(hitpts,[],4)';
    [a,e,~] = getae(flathitpts,1,L,c);
    scatter(a*LinAU,e,'.')
    hold on;
    drawnow;
    
    c = stopCaching(c);
    
end

c = cs(c,'s.i.odeopts',odeset('Events',[],'RelTol',3e-14,'AbsTol',1e-15));

flathitpts = reshape(hitpts,[],4)';
[a,e,omega] = getae(flathitpts,1,L,c);
scatter(a*LinAU,e,'.')