%Attempts to analyze what happens when we make small displacements from the
%singularity.

close all;

%We create a circle of initial conditions that will be used to create the
%positions and velocities.
n = 101;

rdelta = 1e-3;
thetas = linspace(0,2*pi,n);
%The chosen arch intersection energy in the Sun-Jupiter problem
energy = -1.1413;
mass = 2;

clear conditions

for i = 1:n
    conditions(:,i) = energyCorrector(rdelta,thetas(i),energy,mass,c);
    %We multiply the velocities by -1 to get the stable manifold 
    conditions(3:4,i) = conditions(3:4,i)*-1;
end

c = cs(c,'s.i.odeopts',odeset('Events',@xsection,'Reltol',3e-14,...
                              'AbsTol',1e-15));
c = startCaching(c);

endtime = -4.0*pi;

tspan = linspace(0,endtime,1000)

hold on

section = {};

perturbation = 1.1;

for i = 1:size(conditions,2)
    [sol,c] = integ(tspan,conditions(:,i),c);
    
    %We filter out the u4 intersections beyond a certain point
    u4hits = sol.ye(:,sol.ye(1,:) < -2);
    
    for j = 1:size(u4hits,2)
        if j > numel(section)
            section{j} = [];
        end
        section{j} = [u4hits(:,j) section{j}];
    end
    
    %We plot manually; because we needed to get the event crossings, we
    %were unable to use integplot.
    y = deval(sol,tspan);
    p = cplot(y,c,'k');
    p.Color(4) = 0.25;
    drawnow
end
axis equal

%We now create initial conditions that are above and below the manifold

startu = {};
startd = {};
endu = {};
endd = {};
for j = 1:numel(section)
    piece = section{j};
    
    if j > numel(endu)
       startu{j} = [];
       startd{j} = [];
       endu{j} = [];
       endd{j} = [];
    end
    
    meanpt = mean([piece(1,:)
                   piece(3,:)],2)
        
    offsetsect=[piece(1,:)
                piece(3,:)]-meanpt
    
    upoffset = perturbation*offsetsect+meanpt;
    downoffset = (1/perturbation)*offsetsect+meanpt;
            
    for i=1:size(piece,2)

        sectiontraj = piece(:,i);

        up = sectiontraj;
        up(1) = upoffset(1,i);
        up(3) = upoffset(2,i);
        up(5) = 0;

        down = sectiontraj;
        down(1) = downoffset(1,i);
        down(3) = downoffset(2,i);
        down(5) = 0;
        
        startu{j} = [startu{j} up]
        startd{j} = [startd{j} down]

        %We must correct the energies by adjusting vy
        energyHandle = getEnergyHandle(c);

        up(4)=fzero(@(vy)(energyHandle(0,[up(1:3);vy;0])-energy),up(4));
        down(4)=fzero(@(vy)energyHandle(0,[down(1:3);vy;0])-energy,down(4));

        solup=integ(linspace(0,endtime,1000),up,c);
        solun=integ(linspace(0,-1.1*endtime,1000),up,c);

        soldp=integ(linspace(0,endtime,1000),down,c);
        soldn=integ(linspace(0,-1.1*endtime,1000),down,c);

        yup = deval(solup,linspace(0,endtime,1000));
        pup = cplot(yup,c,'b');
        pup.Color(4) = 0.25;

        yun = deval(solun,linspace(0,-1.1*endtime,1000));
        endu{j} = [endu{j} yun(:,end)]
        pun = cplot(yun,c,'b');
        pun.Color(4) = 0.25;

        ydp = deval(soldp,linspace(0,endtime,1000));
        pdp = cplot(ydp,c,'r');
        pdp.Color(4) = 0.25;

        ydn = deval(soldn,linspace(0,-1.1*endtime,1000));
        endd{j} = [endd{j} ydn(:,end)]
        pdn = cplot(ydn,c,'r');
        pdn.Color(4) = 0.25;

        drawnow;

    end
end
%The Sun-Jupiter distance in AU
L = 5.2044;

se1 = figure;
hold on;
se2 = figure;
hold on;
for i = 1:numel(section)
    piece = section{i};
    startupiece = startu{i};
    startdpiece = startd{i};
    endupiece = endu{i};
    enddpiece = endd{i};
    
    figure(se1)
    plot(piece(1,:),piece(3,:),'k.-')
    plot(startupiece(1,:),startupiece(3,:),'b.-')
    plot(startdpiece(1,:),startdpiece(3,:),'r.-')
    
    figure(se2)
    [ap,ep] = getae(piece,1,L,c)
    plot(ap,ep,'k.-')
    [apu,epu] = getae(endupiece,1,L,c)
    plot(apu,epu,'b.')
    [apd,epd] = getae(enddpiece,1,L,c)
    plot(apd,epd,'r.')
    
    
end

% hold on;
% figure(se1)
% scatter(up(1,:),up(3,:),'b')
% scatter(down(1,:),down(3,:),'r')
% 
% figure(se2)
% [au,eu] = getae(up,1,L,c)
% scatter(au,eu,'b')
% [ad,ed] = getae(down,1,L,c)
% scatter(ad,ed,'r')
% %axis equal

c = stopCaching(c);
c = cs(c,'s.i.odeopts',odeset('Events',[],'Reltol',3e-14,'AbsTol',1e-15));