%Attempts to analyze what happens when we make small displacements from the
%singularity.

close all;

%We create a circle of initial conditions that will be used to create the
%positions and velocities.
n = 100;

rdelta = 1e-6;
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

%conditions = conditions(:,50:100);

c = useVelocity(c);

jacobian = getJacobianHandle(cg(c,'d.eqns'),c);

mu = cg(c,'p.mu');
%We look for interceptions within 0.01 and 5 radii of the hill sphere
switchScaling = 0.01;
endScaling = 5;
rh = (mu/3)^(1/3);
switchRadius = switchScaling*rh;
endRadius = endScaling*rh;
primary = 2;
numPts = 10;

%We have to retrieve the f integration equations in order to use
%orbitalelemprop.
eqns = formula(cg(c,'d.eqns'));
feqn = getEquationsHandle(eqns(5),c)

%We send in a custom integrator used for orbital element propagation.
%The hitnum is 1 since we are switching at the first event.
c = cs(c,'d.alt.integ1',@(dummy1,tspan,y0,dummy2)orbitalelemprop(dummy1,...
                          tspan,y0,dummy2,numPts,switchRadius,mu,...
                          primary,1,feqn));


c = cs(c,'s.i.odeopts',odeset('Events',@(t,y)withinRadius_and_xsection(...
       t,y,mu,switchRadius),'Reltol',3e-14,'AbsTol',1e-15));
   
c = startCaching(c);

endtime = -4.0*pi;

tspan = linspace(0,endtime,1000)

hold on

section = {};

%perturbation = 1.00001;
perturbation = 1.00001;

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
    
    %We plot
    y = deval(sol,tspan);
    p = cplot(y,c,'k');
    p.Color(4) = 0.25;
    drawnow
end
axis equal

c = stopCaching(c);

%We create two events: one for switching to the orbital element integrator
%and one for forcing a stop to integration after a certain distance has
%been reached while leaving the vicinity of the second primary.
c = cs(c,'s.i.odeopts',odeset('Events',...
       @(t,y)withinRadius(t,y,cg(c,'p.mu'),[switchRadius;endRadius],...
       [1,0],[-1;1]),'Reltol',3e-14,'AbsTol',1e-15));

c = startCaching(c);

%We now create initial conditions that are above and below the manifold

startu = {};
startd = {};
endu = {};
endd = {};
lengths = {};
for j = 1:1
    piece = section{j};
    
    if j > numel(endu)
       startu{j} = [];
       startd{j} = [];
       endu{j} = [];
       endd{j} = [];
    end
    
    %We calculate the arc length of each point:
    lengthstemp=[];
    lengthstemp(1) = 0;
    for i = 2:size(piece,2)
        lengthstemp(i) = norm(piece(:,i)-piece(:,i-1))+lengthstemp(i-1);
    end
    
    lengths{j} = lengthstemp;
    
    %We now calculate the points inside and outside the manifold.
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
       
        %[pup,c,solup]=integplot(linspace(0,endtime,1000),up,c,'b');
        [pun,c,solun]=integplot(linspace(0,-2*endtime,1000),up,c,'b');

        %[pdp,c,soldp]=integplot(linspace(0,endtime,1000),down,c,'r');
        [pdn,c,soldn]=integplot(linspace(0,-2*endtime,1000),down,c,'r');
        
        if numel(solun)>1
            solunEnd = solun{end};
        else
            solunEnd = solun;
        end
        if ~isempty(solunEnd.ye)
            scatter(solunEnd.ye(1,:),solunEnd.ye(2,:),'b');
            endd{j} = [endd{j} solunEnd.ye(:,end)];
        end
        
        if numel(soldn)>1
            soldnEnd = soldn{end};
        else
            soldnEnd = soldn;
        end
        if ~isempty(soldnEnd.ye)
            scatter(soldnEnd.ye(1,:),soldnEnd.ye(2,:),'r');
            endu{j} = [endu{j} soldnEnd.ye(:,end)];
        end
        
        %pup.Color(4) = 0.25;
        pun.Color(4) = 0.25;
        %pdp.Color(4) = 0.25;
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
for i = 1:1
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
    %[ap,ep] = getae(piece,1,L,c)
    %plot(1:size(piece,2),ep,'k.-')
    [apu,epu] = getae(endupiece,1,L,c)
    plot(lengths{i},epu,'b.')
    [apd,epd] = getae(enddpiece,1,L,c)
    plot(lengths{i},epd,'r.')
    
    
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
c = cs(c,'d.alt.integ1',{})

function [position,isterminal,direction]=withinRadius_and_xsection(t,...
                                            y,mu,switchRadius)
    %Combines withinRadius and xsection in the fashion required above.
    [position1,isterminal1,direction1]=withinRadius(t,y,mu,...
                       switchRadius,1,-1);
    [position2,isterminal2,direction2]=xsection(t,y);
    
    position = [position1; position2];
    isterminal = [isterminal1; isterminal2];
    direction = [direction1; direction2];
    
end