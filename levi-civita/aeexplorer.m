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

aocomega = deg2rad(2.73833e2-360);

eHandle = cg(c,'d.Ehandle');

%The Sun-Jupiter distance in km.
L = 7.784e8;
%The length of an astronomical unit expressed in km.
AUinkm = 149597870.691;
%The Sun-Jupiter distance in AU.
LinAU = L/AUinkm;
%The radius of Jupiter in km.
Jrkm = 71492;
%The radius of Jupiter in denormalized units.
Jr = Jrkm / L

T = 11.829163;
tscale = T/(2*pi);

explorerEndTimeInLC = 10;

mu = cg(c,'p.mu');

%We build out arrays of initial energies and semi-major axes
energyRange = linspace(-1.5194,0,75);
ainAU = linspace(8,20,75);

a = zeros(1,numel(ainAU));
e = zeros(numel(energyRange),numel(ainAU));
mindist = zeros(numel(energyRange),numel(ainAU));
endtimestandard = zeros(numel(energyRange),numel(ainAU));

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
        acur = a(i)
        e(j,i) = fzero(@(eref)(eHandle(0,[aeomega2r(acur,eref,...
                                 aocomega,c); 0])-energy),0.78);
        ic = aeomega2r(a(i),e(j,i),aocomega,c);

        [~,encounteric] = standard2lc(0,ic,mu);

        %We now integrate the trajectory in the Levi-Civita frame
        solbd = integ(linspace(0,explorerEndTimeInLC,1000),encounteric,c);
        [ttrajbd,ytrajbd] = lc2standard(solbd.x,solbd.y,mu);
        
        endtimestandard(j,i) = ttrajbd(end);
        
        %cplot(ytrajbd,c)
        %drawnow;

        %Now, we calculate the minimum position space distance between the
        %trajectory and m_2
        mindist(j,i) = min(vecnorm(ytrajbd(1:2,:) - [1-mu 0]'));
        
    end
    
    c = stopCaching(c);
    
end

flatainAU = reshape((ones(numel(energyRange),1)*ainAU)',1,[]);
flate = reshape(e',1,[]);
flatmindist = reshape(mindist',1,[]);

scatter(flatainAU,flate,100,log(flatmindist),'.');