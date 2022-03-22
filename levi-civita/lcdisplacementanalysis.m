%A displacement analysis that uses integration in the Levi-Civita
%regularization.

close all;

%We create a circle of initial conditions that will be used to create the
%positions and velocities.
n = 20;

rdelta = 1e-8;
%Because of how the LC representation maps back to the standard
%coordinates, we only generate initial conditions between 0 and pi.
thetas = linspace(0,pi,n);
%The chosen arch intersection energy in the Sun-Jupiter problem
energy = -1.1413;


c = cs(c,'p.E',energy);

eHandle = cg(c,'d.Ehandle');

trajfig = figure;
hold on;

endtime = -12;
endRadius = 1;

%We need to square the radius in Levi Civita coords
rh = (cg(c,'p.mu')/3)^(1/3);
scaling = 10;
detectDist = (rh*scaling)^2;

section = {};
wasHit = [];

c = cs(c,'s.i.odeopts',odeset('Events',...
       @(t,y)withinRadius_and_xsection(t,y,detectDist),'RelTol',3e-14,...
                                                       'AbsTol',1e-15));

c = startCaching(c);


for i = 1:n
    %For each initial condition on the manifold, we correct the energy
    %and then integrate in Levi-Civita coordinates
    disp('Condition #:')
    disp(i)
    conditions(:,i) = energyCorrectorLC(rdelta,thetas(i),c);
    %We multiply the velocities by -1 to get the stable manifold 
    conditions(3:4,i) = conditions(3:4,i)*-1;
    [sol,c] = integ(linspace(0, endtime, 2),conditions(:,i),c);
    %[~,~,sol] = integplot(linspace(0, endTime, 1000),conditions(:,i),c);
    
    times = linspace(sol.x(1),sol.x(end),5000)
    
    %We convert the whole trajectory back to standard coordinates:
    [tstrtraj,ystrtraj] = lc2standard(times,deval(times,sol),cg(c,'p.mu'));
    cplot(ystrtraj,c,'k')
    
    %Just as a check, we list the energy after integration. It should be
    %approximately zero.
    disp('Energy at end of trajectory in standard representation:')
    disp(eHandle(tstrtraj(end),[ystrtraj(:,end);tstrtraj(end)]))
    
    %We now convert any hits back to standard coordinates and filter them:
    [hitlct,hitlcy] = lc2standard(sol.xe,sol.ye,cg(c,'p.mu'));
    
    %These are the outgoing hits for determining eccentricity:
    hitlcyoutgoing(:,i) = hitlcy(:,sol.ie == 1);
    
    %These are the hits on the Poincare section:
    hitlcyfiltered = hitlcy(:,(hitlcy(1,:) < -2) & (sol.ie == 2));
    
    wasHit(i) = ~isempty(hitlcyfiltered);
    
    cplot(hitlcyfiltered,c,'ok');
    
    %We assemble the section points
    for j = 1:size(hitlcyfiltered,2)
        if j > numel(section)
            section{j} = [];
        end
        section{j} = [hitlcyfiltered(:,j) section{j}];
        
    end
    
    disp('------------')
   
    drawnow;
end

c = stopCaching(c);   

c = cs(c,'s.i.odeopts',odeset('Events',@(t,y)withinRadius(t,y,...
                       0,detectDist,1,1),'RelTol',3e-14,...
                       'AbsTol',1e-15));

c = startCaching(c);

%We now create initial conditions that are above and below the manifold

startu = {};
startd = {};
endu = {};
endd = {};
lengths = {};

perturbation = 1e-5;

%We now create initial conditions that are above and below the manifold
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
    
    proj = [piece(1,:)
     piece(3,:)];
    
    %We now calculate the points inside and outside the manifold.
    meanpt = mean(proj,2)
        
    %We calculate and then normalize the vector pointing from the center of
    % the manifold to each point on it. 
    offsetsect = proj-meanpt;
    offsetsect = offsetsect./vecnorm(offsetsect);  
    
%     numManifoldPts = size(proj,2);
%     
%     offsetsect = [];
%     
%     for i = 1:numManifoldPts
%         %We approximate the tangent vector using the points surrounding the
%         %point under consideration
%         if i > 1
%             prevpt = proj(:,i-1);
%         else
%             prevpt = proj(:,end);
%         end
%         
%         if i < numManifoldPts
%             nextpt = proj(:,i+1);
%         else
%             nextpt = proj(:,1);
%         end
%         
%         tangentVec = nextpt - prevpt;
%         perpendicularVec = [-tangentVec(2)
%                              tangentVec(1)];
%         offsetsect = [offsetsect perpendicularVec];
%         
%     end    
            
    upoffset = proj + perturbation.*offsetsect;
    downoffset = proj - perturbation.*offsetsect;
            
    for i=1:size(piece,2)

        sectiontraj = piece(:,i);

        up = sectiontraj;
        up(1) = upoffset(1,i);
        up(3) = upoffset(2,i);


        down = sectiontraj;
        down(1) = downoffset(1,i);
        down(3) = downoffset(2,i);

        
        startu{j} = [startu{j} up]
        startd{j} = [startd{j} down]

        up(4)=fzero(@(py)(eHandle(0,[up(1:3);py;0])-energy),up(4));
        down(4)=fzero(@(py)(eHandle(0,[down(1:3);py;0])-energy),down(4));
        
        %Now, we must convert back to Levi-Civita form:
        [~,uplc] = standard2lc(zeros(1,size(up,2)),up,cg(c,'p.mu'));
        [~,downlc] = standard2lc(zeros(1,size(down,2)),down,cg(c,'p.mu'));
       
        %[pup,c,solup]=integplot(linspace(0,endtime,1000),up,c,'b');
        solun=integ(linspace(0,-5*endtime,1000),uplc,c);

        %[pdp,c,soldp]=integplot(linspace(0,endtime,1000),down,c,'r');
        soldn=integ(linspace(0,-5*endtime,1000),downlc,c);
        
        if numel(solun)>1
            solunEnd = solun{end};
        else
            solunEnd = solun;
        end
        
        timesun = linspace(solunEnd.x(1),solunEnd.x(end),50000);
        [~,trajunpts] = lc2standard(timesun,deval(timesun,solunEnd),...
                                            cg(c,'p.mu'));
        cplot(trajunpts,c,'b')
        
        if ~isempty(solunEnd.ye)
            [endTimes,unendPts] = lc2standard(solunEnd.xe,solunEnd.ye,...
                                            cg(c,'p.mu'));
            
            cplot(unendPts,c,'bo');
            endu{j} = [endu{j} unendPts(:,end)];
        end
        
        if numel(soldn)>1
            soldnEnd = soldn{end};
        else
            soldnEnd = soldn;
        end
        
        timesdn = linspace(soldnEnd.x(1),soldnEnd.x(end),50000);
        [~,trajdnpts] = lc2standard(timesdn,deval(timesdn,soldnEnd),...
                                            cg(c,'p.mu'));
        cplot(trajdnpts,c,'r')
        
        if ~isempty(soldnEnd.ye)
            [endTimes,dnendPts] = lc2standard(soldnEnd.xe,soldnEnd.ye,...
                                            cg(c,'p.mu'));
            
            cplot(dnendPts,c,'ro');
            endd{j} = [endd{j} dnendPts(:,end)];
            
        end
        
%         %pup.Color(4) = 0.25;
%         pun.Color(4) = 0.25;
%         %pdp.Color(4) = 0.25;
%         pdn.Color(4) = 0.25;

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
    [ap,ep] = getae(hitlcyoutgoing,1,L,c)
    plot(lengths{i},ep(logical(fliplr(wasHit))),'k.-')
    %plot(1:numel(ep),ep,'k.-')
    [apu,epu] = getae(endupiece,1,L,c)
    plot(lengths{i},epu,'b.')
    [apd,epd] = getae(enddpiece,1,L,c)
    plot(lengths{i},epd,'r.')
    
    
end

c = stopCaching(c);   

c = cs(c,'s.i.odeopts',odeset('Events',[],'RelTol',3e-14,'AbsTol',1e-15));

function [position,isterminal,direction]=withinRadius_and_xsection(t,...
                                            y,radDist)
    %Combines withinRadius and xsection in the fashion required above.
    [position1,isterminal1,direction1]=withinRadius(t,y,0,...
                       radDist,0,1);
    [position2,isterminal2,direction2]=xsection(t,y);
    
    position = [position1; position2];
    isterminal = [isterminal1; isterminal2];
    direction = [direction1; direction2];
    
end