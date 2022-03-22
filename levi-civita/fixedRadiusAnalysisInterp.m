%An analysis of initial conditions on different sides of the manifold
%at a fixed radius. Uses interpolated initial conditions for things to get
%more precise and less noisy results.

close all;

%The chosen arch intersection energy in the Sun-Jupiter problem
energy = -1.1413;
c = cs(c,'p.E',energy);
eHandle = cg(c,'d.Ehandle');

%The Sun-Jupiter distance in km.
L = 7.784e8;
%The radius of Jupiter in km.
Jrkm = 71492;
%The radius of Jupiter in denormalized units.
Jr = Jrkm / L

mu = cg(c,'p.mu');

%We create 2n initial conditions; n are on one side of the manifold, and n
%are on the other side of it.
n = 1000;

%This is the radius at which the initial conditions will be generated; it
%corresponds to the closest encounter with the singularity.
%rpd = Jr/20;
%rpd = 0.5030690810017047*Jr;
rpd = 0.5030690810017048*Jr;

%We generate an array of different theta values, corresponding to different
%points of tangency on the periapse circle. We generate an extra value at 
%2*pi and then remove it because it is identical to the theta = 0 case,
%leaving us with n unique values (mod 2*pi of course).
theta = linspace(0,2*pi,n+1);
theta = theta(1:end-1);

%mindist
%theta = theta(139);
%maxdist
%theta = theta(386);

%theta = theta([139,250,386]);

%We continue the convention of calling the two sets of initial conditions
%"up" and "down."
up = zeros(4,n);
down = zeros(4,n);

%This is the maximum amount of time to integrate, but we expect that we
%should not have to integrate for this long unless something goes wrong (as
%all of these initial conditions should intercept the detection surface of
%section).
endtime = -12;

%We construct the detection surface of section. We need to square the
%desired radius (in standard coordinates) to express it properly in 
%Levi-Civita coords
rh = (mu/3)^(1/3);
scaling = 12%1.44%15;
detectDist = (rh*scaling)^2;

%c = cs(c,'s.o.v.dmode','velocity');

c = cs(c,'s.i.odeopts',odeset('Events', @(t,y)withinRadius(t,y,0,...
                       detectDist,1,1),'RelTol',3e-14,'AbsTol',1e-15));

c = startCaching(c);

%We call getFromCache for the sole purpose of caching the equations of
%motion.
[~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

uplclist = zeros(4,numel(theta));
downlclist = zeros(4,numel(theta));
endufd = zeros(4,numel(theta));
endubd = zeros(4,numel(theta));
enddbd = zeros(4,numel(theta));

soldnfdEnd = {};
soldnbdEnd = {};
solunfdEnd = {};
solunbdEnd = {};

parfor i = 1:numel(theta)
    
    disp(i)
    
    %For each initial condition, we create a function handle to be passed
    %to fzero in order to refine a suitable thetadot for the chosen energy,
    %and then we pass it. 
    refineHandle = @(thetadot)(eHandle(0,[m2polar2cart(...
                                          [rpd
                                           theta(i)
                                           0
                                           thetadot],mu);0])-energy);
    %We take the absolute value for consistency.
    thetadot = abs(fzero(refineHandle,1));
    
    up(:,i) = m2polar2cart([rpd theta(i) 0 thetadot].',mu);
    down(:,i) = m2polar2cart([rpd theta(i) 0 -thetadot].',mu);
    
    calcSepAngle(rpd,thetadot,mu)
    
    %Now, we must convert back to Levi-Civita form:
    [~,uplc] = standard2lc(0,up(:,i),mu);
    [~,downlc] = standard2lc(0,down(:,i),mu);
    
    uplclist(:,i) = up(:,i);
    downlclist(:,i) = down(:,i);

    solunfd=integ(linspace(0,5*endtime,1000),uplc,c);
    soldnfd=integ(linspace(0,5*endtime,1000),downlc,c);
    solunbd=integ(linspace(0,-5*endtime,1000),uplc,c);
    soldnbd=integ(linspace(0,-5*endtime,1000),downlc,c);

    %%%% Handling the results of integrating up points forward
    if numel(solunfd)>1
        solunfdEnd{i} = solunfd{end};
    else
        solunfdEnd{i} = solunfd;
    end

    if ~isempty(solunfdEnd{i}.ye)
        [endTimes,unfdendPts] = lc2standard(solunfdEnd{i}.xe,solunfdEnd{i}.ye,...
                                        mu);

        endufd(:,i) = unfdendPts(:,end);
    end

    %%%% Handling the results of integrating up points backward
    if numel(solunbd)>1
        solunbdEnd{i} = solunbd{end};
    else
        solunbdEnd{i} = solunbd;
    end

    if ~isempty(solunbdEnd{i}.ye)
        [endTimes,unbdendPts] = lc2standard(solunbdEnd{i}.xe,solunbdEnd{i}.ye,...
                                        mu);

        endubd(:,i) = unbdendPts(:,end);
    end

    %%%% Handling the results of integrating down points backward
    if numel(soldnbd)>1
        soldnbdEnd{i} = soldnbd{end};
    else
        soldnbdEnd{i} = soldnbd;
    end

    if ~isempty(soldnbdEnd{i}.ye)
        [endTimes,dnbdendPts] = lc2standard(soldnbdEnd{i}.xe,soldnbdEnd{i}.ye,...
                                        mu);
        enddbd(:,i) = dnbdendPts(:,end);

    end

end

%Now, we need to sort the points according to their angles along
%the detection circle.
[endubdsorted, endubdi, endubdthetas] = angleSort(endubd,mu);
endufdsorted = endufd(:,endubdi);
endufdpolar = cart2m2polar(endufdsorted, mu);
endufdthetas = endufdpolar(2,:);

[enddbdsorted, enddbdi, enddbdthetas] = angleSort(enddbd,mu);

%We interpolate the down points integrated backwards in order to find down
%points which more exactly correspond with up points.
enddbdint = spline(enddbdthetas,enddbdsorted,endubdthetas);
enddbdintpolar = cart2m2polar(enddbdint, mu);
enddbdintthetas = enddbdintpolar(2,:);

intdownlclist = zeros(4,numel(theta));
enddfdint = zeros(4,numel(theta));
solintdnfdEnd = {};

%We now integrate the interpolated down intercepts back through the origin
%and to the detection circle once more.
parfor i = 1:numel(enddbdintthetas)
    
    %Now, we must convert back to Levi-Civita form:
    [~,intdownlc] = standard2lc(0,enddbdint(:,i),mu);
    
    intdownlclist(:,i) = intdownlc;

    solintdnfd=integ(linspace(0,10*endtime,1000),intdownlc,c);

    %%%% Handling the results of integrating down points forward
    if numel(solintdnfd)>1
        solintdnfdEnd{i} = solintdnfd{end};
    else
        solintdnfdEnd{i} = solintdnfd;
    end

    if ~isempty(solintdnfdEnd{i}.ye)
        [endTimes,dnfdendPts] = lc2standard(solintdnfdEnd{i}.xe,solintdnfdEnd{i}.ye,...
                                        mu);

        enddfdint(:,i) = dnfdendPts(:,end);
    end
end

c = stopCaching(c);

%We get the theta values:
enddfdintpolar = cart2m2polar(enddfdint, mu);
enddfdintthetas = enddfdintpolar(2,:);

%We get the predist and the postdist:
predist = vecnorm(endubdsorted - enddbdint).';
postdist = vecnorm(endufdsorted - enddfdint).';

figure;
hold on;

% [apu,epu] = getae(uplclist,1,L,c);
% plot(theta,epu,'b.-')
% [apd,epd] = getae(downlclist,1,L,c);
% plot(theta,epd,'r.-')

%We plot the eccentricity plots:
[apufd,epufd] = getae(endufdsorted,1,L,c);
plot(endubdthetas,epufd,'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Post-encounter eccentricity (+ trajectories)')

figure;
hold on;

[apufd,epufd] = getae(endubdsorted,1,L,c);
plot(endubdthetas,epufd,'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Pre-encounter eccentricity (+ trajectories)')

figure;
hold on;

[apdfd,epdfd] = getae(enddfdint,1,L,c);
plot(enddbdintthetas,epdfd,'r-')
xlabel('Pre-encounter theta (Interpolated - trajectories)')
ylabel('Post-encounter eccentricity (Interpolated - trajectories)')

figure;
hold on;

[apdfd,epdfd] = getae(enddbdsorted,1,L,c);
plot(enddbdintthetas,epdfd,'r-')
xlabel('Pre-encounter theta (Interpolated - trajectories)')
ylabel('Pre-encounter eccentricity (Interpolated - trajectories)')

%We also plot the delta-t plots:
figure;
hold on;
plot(endubdthetas,mod(endufdthetas-endubdthetas,2*pi),'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Change in theta after encounter (+ trajectories)')

figure;
hold on;
plot(enddbdthetas,enddfdintthetas-enddbdintthetas,'b-')
xlabel('Pre-encounter theta (- trajectories)')
ylabel('Change in theta after encounter (- trajectories)')

%We verify the matched points:
%figure;
%hold on;
%cplot(endubdsorted(:,1:floor(end/2)),c,'o')
%cplot(enddbdsorted(:,matchindex(1:floor(end/2))),c,'o')

%We plot the magnitude of the separation before and after the encounter
%against the integrated backwards theta value of the up point:
figure;
hold on;
plot(endubdthetas, predist)
plot(endubdthetas, postdist)
legend('Pre-encounter dist','Post-encounter dist')

%We also plot the magnitude of the multiplier k, which is defined as
%postdist / predist:
figure; 
plot(endubdthetas,postdist./predist)

c = cs(c,'s.i.odeopts',odeset('Events', [],'RelTol',3e-14,'AbsTol',1e-15));