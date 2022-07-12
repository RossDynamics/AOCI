%Attempts to replicate the "arches of chaos" plot by plotting close
%encounter distance. SUN_JUPITER_ER3BP_CONTEXT_INTEGTIME should be used to
%create the context object.
%ic = aeomega2r(11.2/LinAU,0.7845,aocomega,c)
%ic = aeomega2r(11.21/LinAU,0.7846,aocomega,c)
%%%%%%%%ic = aeomega2r(11.213/LinAU,0.78465,aocomega,c)
%ic = aeomega2r(11.215/LinAU,0.7847,aocomega,c)
%ic = aeomega2r(11.22/LinAU,0.78475,aocomega,c)
%ic = aeomega2r(11.3/LinAU,0.7856,aocomega,c)

%ic = aeomega2r(12/LinAU,0.79295,aocomega,c)
%ic = aeomega2r(12.695/LinAU,0.8,aocomega,c)
%ic = aeomega2r(12.59/LinAU,0.79895,aocomega,c)
%ic = aeomega2r(12.8/LinAU,0.80105,aocomega,c)
%ic = aeomega2r(13/LinAU,0.80305,aocomega,c)
%ic = aeomega2r(13.3/LinAU,0.80595,aocomega,c)
%ic = aeomega2r(13.4/LinAU,0.8069,aocomega,c)
%ic = aeomega2r(13.48/LinAU,0.80765,aocomega,c)
%%%%%%%%ic = aeomega2r(13.49/LinAU,0.80775,aocomega,c)
%ic = aeomega2r(13.5/LinAU,0.80785,aocomega,c)
%ic = aeomega2r(14/LinAU,0.8125,aocomega,c)

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
tscale = T/(2*pi);

%Our epoch is 30 September 2012
epochinyears = years(datetime(2012,09,30)-datetime(2000,1,1,12,0,0));
epochdenorm = epochinyears/tscale;

%The J2000 mean longitude angle:
j2000Angle = deg2rad(34.40438);

%epochAngle = epochdenorm + j2000Angle;
epochAngle = 0.8159069576305198;

%This is used to put an upper bound on how long the integrator is allowed
%to integrate before it is interrupted. It should ideally be rather large. 
explorerEndTimeInLC = 5000;
%This is the actual end time in nondimensional time units.
standardEndTime = 78 / tscale;

mu = cg(c,'p.mu');

%We build out arrays of initial energies and semi-major axes
energyRange = linspace(-1.5194,0,200);
%ainAU = linspace(8,20,100);
ainAU = linspace(8,20,400);

a = zeros(1,numel(ainAU));
e = zeros(numel(energyRange),numel(ainAU));
mindist = zeros(numel(energyRange),numel(ainAU));
endtimestandard = zeros(numel(energyRange),numel(ainAU));

c = cs(c,'s.i.odeopts',odeset('Events',@(t,y)untilStandardTime(t,y,...
              standardEndTime),'RelTol',3e-14,'AbsTol',1e-15));

for j = 1:numel(energyRange)

    energy = energyRange(j);

    c = cs(c,'p.E',energy);

    c = startCaching(c);

    %We call getFromCache for the sole purpose of caching the equations of
    %motion.
    [~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

    disp(energy)
    
    parfor i = 1:numel(ainAU)
        %disp(j)
        %disp(energy)
        %disp(i)
        %disp(ainAU(i))

        %We need to refine the eccentricity that fixes the energy for the
        %chosen semi-major axis. Semi-major axis doesn't vary between
        %energies, but eccentricity does, and so the array dimensions
        %reflect this.
        a(i) = ainAU(i)/LinAU;
        acur = a(i);
        %e(j,i) = fzero(@(eref)(eHandle(0,[aeomega2r(acur,eref,...
        %                         aocomega,c); 0])-energy),0.78);
        %ic = aeomega2r(a(i),e(j,i),aocomega,c);
        e(j,i) = fzero(@(eref)(eHandle(0,[aeomegaM2r(acur,eref,...
                                 aocomega,M,epochAngle,c); 0])-energy),0.78);
        ic = aeomegaM2r(a(i),e(j,i),aocomega,M,epochAngle,c);
        
        [~,~,omega] = getae(ic,1,1,c);
        disp('Omega:')
        disp(omega)
        
        [~,encounteric] = standard2lc(0,ic,mu);

        %We now integrate the trajectory in the Levi-Civita frame
        try
            solbd = integ(linspace(0,explorerEndTimeInLC,1000),[encounteric;0],c);
        catch exception
            disp('An exception occurred while integrating. Ignoring.')
            disp('Exception details:')
            exception
            continue;
        end
        [ttrajbd,ytrajbd] = lc2standard(solbd.x,solbd.y,mu);
        
        %We permit NaN trajectories to pass silently
        if ~numel(solbd.ie) && ~any(isnan(solbd.y),'all')
            throw(MException('AOCI:ExplorerEndTimeInLCTooSmall',...
'explorerEndTimeInLC needs to be larger so that integration can finish.'));
        end
        
        endtimestandard(j,i) = ttrajbd(end);
        
        %cplot(ytrajbd,c)
        %drawnow;

        %Now, we calculate the minimum position space distance between the
        %trajectory and m_2
        mindist(j,i) = min(vecnorm(ytrajbd(1:2,:) - [1-mu 0]'));
        
    end
    
    c = stopCaching(c);
    
end

c = cs(c,'s.i.odeopts',odeset('Events',[],'RelTol',3e-14,'AbsTol',1e-15));

flatainAU = reshape((ones(numel(energyRange),1)*ainAU)',1,[]);
flate = reshape(e',1,[]);
flatmindist = reshape(mindist',1,[]);

imagesc([2 20], [0,1],flipud(aoc78yrs))
hold on;
scatter(flatainAU,flate,100,log(1./flatmindist),'.');
set(gca,'ydir','normal');
colormap('bone')

%scatter(flatainAU,flate,100,log(flatmindist),'.');