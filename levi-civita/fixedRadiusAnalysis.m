%An analysis of initial conditions on different sides of the manifold
%at a fixed radius. If this script is being called by
%calcPatchedFullError, then neither rpd and scaling should be initialized!

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
n = 400;

%We construct the detection surface of section. We need to take the 
%square root of the desired radius (in standard coordinates) to express it
%properly in Levi-Civita coords
rh = (mu/3)^(1/3);
%scaling = 7;
%scaling = 1;%8%1.44%15;
detectDistStd = rh*scaling;
detectDist = sqrt(detectDistStd);

%This is the radius at which the initial conditions will be generated; it
%corresponds to the closest encounter with the singularity.
%rpd = Jr;
%rpd = 0.5030690810017047*Jr;
%rpd = 0.5030690810017048*Jr;
%rpd =Jr/10;
%rpd = Jr/40;
%rpd = Jr/180;


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

%c = cs(c,'s.o.v.dmode','velocity');

c = cs(c,'s.i.odeopts',odeset('Events', @(t,y)withinRadius(t,y,0,...
                       detectDist,1,1),'RelTol',3e-14,'AbsTol',1e-15));

c = startCaching(c);

%We call getFromCache for the sole purpose of caching the equations of
%motion.
[~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

plusthetadotlist = zeros(1,numel(theta));
minusthetadotlist = zeros(1,numel(theta));

uplclist = zeros(4,numel(theta));
downlclist = zeros(4,numel(theta));
endufd = zeros(4,numel(theta));
enddfd = zeros(4,numel(theta));
endubd = zeros(4,numel(theta));
enddbd = zeros(4,numel(theta));

soldnfdEnd = {};
soldnbdEnd = {};
solunfdEnd = {};
solunbdEnd = {};

%We initialize the guess
plusguess = 1e+5;
minusguess = -1e+5;

for i = 1:numel(theta)
    %For each initial condition, we create a function handle to be passed
    %to fzero in order to refine a suitable thetadot for the chosen energy,
    %and then we pass it. 
    refineHandle = @(thetadot)(eHandle(0,[m2polar2cart(...
                                          [rpd
                                           theta(i)
                                           0
                                           thetadot],mu);0])-energy);
    
                                       
    %We take the absolute value for consistency.
    plusthetadot = fzero(refineHandle,plusguess);
    minusthetadot = fzero(refineHandle,minusguess);
    
    plusthetadotlist(i) = plusthetadot;
    minusthetadotlist(i) = minusthetadot;
    
    %We this thetadot as a guess for the next one.
    plusguess = plusthetadot;
    minusguess = minusthetadot;
        
    up(:,i) = m2polar2cart([rpd theta(i) 0 plusthetadot].',mu);
    down(:,i) = m2polar2cart([rpd theta(i) 0 minusthetadot].',mu);
    
    %calcSepAngle(rpd,thetadot,mu)
    
    %Now, we must convert back to Levi-Civita form:
    [~,uplc] = standard2lc(0,up(:,i),mu);
    [~,downlc] = standard2lc(0,down(:,i),mu);
    
    uplclist(:,i) = uplc;
    downlclist(:,i) = downlc;

end

disp('Finished building initial conditions.')

parfor i = 1:numel(theta)
    
    disp(i)
    
    solunfd=integ(linspace(0,5*endtime,1000),uplclist(:,i),c);
    soldnfd=integ(linspace(0,5*endtime,1000),downlclist(:,i),c);
    solunbd=integ(linspace(0,-5*endtime,1000),uplclist(:,i),c);
    soldnbd=integ(linspace(0,-5*endtime,1000),downlclist(:,i),c);

    %%%% Handling the results of integrating up points forward
    if numel(solunfd)>1
        solunfdEnd{i} = solunfd{end};
    else
        solunfdEnd{i} = solunfd;
    end

    timesunfd = linspace(solunfdEnd{i}.x(1),solunfdEnd{i}.x(end),50000);
    [~,trajunfdpts] = lc2standard(timesunfd,deval(timesunfd,solunfdEnd{i}),...
                                        mu);
    %cplot(trajunfdpts,c,'b')

    if ~isempty(solunfdEnd{i}.ye)
        [endTimes,unfdendPts] = lc2standard(solunfdEnd{i}.xe,solunfdEnd{i}.ye,...
                                        mu);

        %cplot(unfdendPts,c,'bo');
        endufd(:,i) = unfdendPts(:,end);
    end

    %%%% Handling the results of integrating down points forward
    if numel(soldnfd)>1
        soldnfdEnd{i} = soldnfd{end};
    else
        soldnfdEnd{i} = soldnfd;
    end

    timesdnfd = linspace(soldnfdEnd{i}.x(1),soldnfdEnd{i}.x(end),50000);
    [~,trajdnfdpts] = lc2standard(timesdnfd,deval(timesdnfd,soldnfdEnd{i}),...
                                        mu);
    %cplot(trajdnfdpts,c,'r')

    if ~isempty(soldnfdEnd{i}.ye)
        [endTimes,dnfdendPts] = lc2standard(soldnfdEnd{i}.xe,soldnfdEnd{i}.ye,...
                                        mu);

        %cplot(dnfdendPts,c,'ro');
        enddfd(:,i) = dnfdendPts(:,end);

    end
    
    %drawnow;
    
    %%%% Handling the results of integrating up points backward
    if numel(solunbd)>1
        solunbdEnd{i} = solunbd{end};
    else
        solunbdEnd{i} = solunbd;
    end

    timesunbd = linspace(solunbdEnd{i}.x(1),solunbdEnd{i}.x(end),1000);
    [~,trajunbdpts] = lc2standard(timesunbd,deval(timesunbd,solunbdEnd{i}),...
                                        mu);
    %cplot(trajunbdpts,c,'b')

    if ~isempty(solunbdEnd{i}.ye)
        [endTimes,unbdendPts] = lc2standard(solunbdEnd{i}.xe,solunbdEnd{i}.ye,...
                                        mu);

        %cplot(unbdendPts,c,'bo');
        endubd(:,i) = unbdendPts(:,end);
    end

    %%%% Handling the results of integrating down points backward
    if numel(soldnbd)>1
        soldnbdEnd{i} = soldnbd{end};
    else
        soldnbdEnd{i} = soldnbd;
    end

    timesdnbd = linspace(soldnbdEnd{i}.x(1),soldnbdEnd{i}.x(end),1000);
    [~,trajdnbdpts] = lc2standard(timesdnbd,deval(timesdnbd,soldnbdEnd{i}),...
                                        mu);
    %cplot(trajdnbdpts,c,'r')

    if ~isempty(soldnbdEnd{i}.ye)
        [endTimes,dnbdendPts] = lc2standard(soldnbdEnd{i}.xe,soldnbdEnd{i}.ye,...
                                        mu);

        %cplot(dnbdendPts,c,'ro');
        enddbd(:,i) = dnbdendPts(:,end);

    end

%         %pup.Color(4) = 0.25;
%         pun.Color(4) = 0.25;
%         %pdp.Color(4) = 0.25;
%         pdn.Color(4) = 0.25;

    %drawnow;
    axis equal;
end

c = stopCaching(c);

%Now, we need to sort the points according to their angles along
%the detection circle.
[endubdsorted, endubdi, endubdthetas] = angleSort(endubd,mu);
endufdsorted = endufd(:,endubdi);
endufdpolar = cart2m2polar(endufdsorted, mu);
endufdthetas = endufdpolar(2,:);

[enddbdsorted, enddbdi, enddbdthetas] = angleSort(enddbd,mu);
enddfdsorted = enddfd(:,enddbdi); 
enddfdpolar = cart2m2polar(enddfdsorted, mu);
enddfdthetas = enddfdpolar(2,:);
%[~, ~, enddfdthetas] = angleSort(enddfdsorted,mu);

%We now assign matched points by ascertaining which down points integrated
%backwards are closest to which up points integrated backwards. predist
%is the initial separation between these points before the encounter, and
%postdist is the separation after the encounter.
predist = [];
postdist = [];
matchindex = [];

enddbdsorted(:,:)

for i = 1:numel(theta)
    searchthrough = enddbdsorted(:,:);
%     if i ~= numel(theta)
%     searchthrough = enddbdsorted(:,i+1:end);
%     else
%         searchthrough = enddbdsorted(:,end);
%     end
    query = endubdsorted(:,i);
    [matchindexi, predisti] = dsearchn(searchthrough.',query.');
    matchindex(i) = matchindexi;
    predist(i) = predisti;
    postdist(i) = vecnorm(endufdsorted(:,i) - enddfdsorted(:,matchindexi)).';
end

figure;
hold on;

% [apu,epu] = getae(uplclist,1,L,c);
% plot(theta,epu,'b.-')
% [apd,epd] = getae(downlclist,1,L,c);
% plot(theta,epd,'r.-')

%We plot the eccentricity plots:
[apufd,epufd,omegaufd,energiesufd] = getae(endufdsorted,1,L,c);
plot(endubdthetas,epufd,'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Post-encounter eccentricity (+ trajectories)')

figure;
hold on;

[apubd,epubd,omegaubd,energiesubd] = getae(endubdsorted,1,L,c);
plot(endubdthetas,epubd,'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Pre-encounter eccentricity (+ trajectories)')

figure;
hold on;

[apdfd,epdfd,omegadfd,energiesdfd] = getae(enddfdsorted,1,L,c);
plot(enddbdthetas,epdfd,'r-')
xlabel('Pre-encounter theta (- trajectories)')
ylabel('Post-encounter eccentricity (- trajectories)')

figure;
hold on;

[apdbd,epdbd,omegadbd,energiesdbd] = getae(enddbdsorted,1,L,c);
plot(enddbdthetas,epdbd,'r-')
xlabel('Pre-encounter theta (- trajectories)')
ylabel('Pre-encounter eccentricity (- trajectories)')

%We also plot the delta-t plots:
figure;
hold on;
plot(endubdthetas,mod(endufdthetas-endubdthetas,2*pi),'b-')
xlabel('Pre-encounter theta (+ trajectories)')
ylabel('Change in theta after encounter (+ trajectories)')

figure;
hold on;
plot(enddbdthetas,mod(enddfdthetas-enddbdthetas,2*pi),'b-')
xlabel('Pre-encounter theta (- trajectories)')
ylabel('Change in theta after encounter (- trajectories)')

%We verify the matched points:
figure;
hold on;
cplot(endubdsorted(:,1:floor(end/2)),c,'o')
cplot(enddbdsorted(:,matchindex(1:floor(end/2))),c,'o')

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
plot(enddbdthetas,postdist./predist)

%We plot comparisons of the semi-major axes before the encounter and
%after the encounter
figure;
hold on;
plot(endubdthetas,apubd)
plot(endubdthetas,apufd)
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('+ trajectories'' SMA w.r.t primary before encounter',...
       '+ trajectories'' SMA w.r.t primary after encounter')

figure;
hold on;
plot(enddbdthetas,apdbd)
plot(enddbdthetas,apdfd)
xlabel('$\theta^-_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('- trajectories'' SMA w.r.t primary before encounter',...
       '- trajectories'' SMA w.r.t primary after encounter')
   
%We plot comparisons of the argument of perigee before the encounter and
%after the encounter
figure;
hold on;
plot(endubdthetas,omegaubd)
plot(endubdthetas,omegaufd)
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('+ trajectories'' AOP w.r.t primary before encounter',...
       '+ trajectories'' AOP w.r.t primary after encounter')

figure;
hold on;
plot(enddbdthetas,omegadbd)
plot(enddbdthetas,omegadfd)
xlabel('$\theta^-_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('- trajectories'' AOP w.r.t primary before encounter',...
       '- trajectories'' AOP w.r.t primary after encounter')
   
%We plot the "spikes" in the semi-major axis graph
[~,examplebd] = max(apubd);
%examplebd = 300;
[~,examplelc] = standard2lc(0,endubdsorted(:,examplebd),mu);

examplesol=integ(linspace(0,6*endtime,2000),examplelc,c);
timesexample = linspace(examplesol.x(1),examplesol.x(end),2000);
[exampletimesstandard,examplepts] = lc2standard(timesexample,deval(timesexample,examplesol),mu);
                                
exampleplot = figure;
hold on;
cplot(examplepts,c);

figure;
[apuexample,epuexample,omegauexample] = getae(examplepts,1,L,c);
plot(timesexample,apuexample)
hold on;
[examplemaxval,examplemaxi] = max(apuexample);
scatter(timesexample(examplemaxi),examplemaxval)

figure(exampleplot);
cplot(examplepts(:,examplemaxi),c,'o')

%We plot comparisons of the semi-major axes before the encounter and
%after the encounter
figure;
hold on;
plot(endubdthetas,apubd)
plot(endubdthetas,apufd)
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('+ trajectories'' SMA w.r.t primary before encounter',...
       '+ trajectories'' SMA w.r.t primary after encounter')

figure;
hold on;
plot(enddbdthetas,apdbd)
plot(enddbdthetas,apdfd)
xlabel('$\theta^-_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('- trajectories'' SMA w.r.t primary before encounter',...
       '- trajectories'' SMA w.r.t primary after encounter')
   
%%

c = cs(c,'s.i.odeopts',odeset('Events', @(t,y)withinRadius(t,y,0,...
                       4*detectDist,1,1),'RelTol',3e-14,'AbsTol',1e-15));

%We plot comparisons of the Keplerian energy before the encounter and
%after the encounter
[~,highlighti]=min(abs(energiesubd-energiesufd))
%highlighti = 170;
%highlighti = 50;
figure;
subplot(2,2,1);
hold on;
plot(endubdthetas,energiesubd,'b.')
plot(endubdthetas,energiesufd,'b-')
%plot(enddbdthetas,energiesdbd,'r.')
plot(enddbdthetas,energiesdfd,'r-')
scatter(endubdthetas(highlighti),energiesubd(highlighti),'bo')
scatter(endubdthetas(highlighti),energiesdbd(highlighti),'ro')
scatter(endubdthetas(highlighti),energiesufd(highlighti),'bo')
scatter(endubdthetas(highlighti),energiesdfd(highlighti),'ro')
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('$E^\pm_{pre}$ w.r.t $m_1$',...
       '$E^+_{post}$ w.r.t $m_1$',... %'$E^-_{pre}$ w.r.t $m_1$',...
       '$E^-_{post}$ w.r.t $m_1$',...
       'Highlighted + trajectory','Highlighted - trajectory',...
       'Interpreter','Latex')
fplot(@(x)0,'k-')
title('(a)','Interpreter','Latex')
%set(gca,'FontSize',20);

%We plot the highlighted trajectory in the rotating frame
subplot(2,2,2);
[~,highlightlcu] = standard2lc(0,endubdsorted(:,highlighti),mu);
[~,highlightlcd] = standard2lc(0,enddbdsorted(:,highlighti),mu);

highlightsolu=integ(linspace(0,3*endtime,2000),highlightlcu,c);
highlightsold=integ(linspace(0,3*endtime,2000),highlightlcd,c);
timeshighlightu = linspace(highlightsolu.x(1),highlightsolu.x(end),2000);
timeshighlightd = linspace(highlightsold.x(1),highlightsold.x(end),2000);
[highlightstdtimesu,highlightptsu] = lc2standard(timeshighlightu,...
                                  deval(timeshighlightu,highlightsolu),mu);
[highlightstdtimesd,highlightptsd] = lc2standard(timeshighlightd,...
                                  deval(timeshighlightd,highlightsold),mu);
cplot(highlightptsu,c,'b-')
hold on;
cplot(highlightptsd,c,'r-')
cplot(highlightptsu(:,1),c,'ko')
cplot(highlightptsd(:,1),c,'ko')
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
legend('Highlighted + trajectory',...
       'Highlighted - trajectory',...
       'Backwards initial condition','Interpreter','Latex')
axis square;
%set(gca,'FontSize',20);
set(gca,'YAxisLocation','right');
title('(b)','Interpreter','Latex')

%We plot the evolution of the argument of perigee
subplot(2,2,3);
hold on;
[~,discontubd] = max(abs(diff(omegaubd)));
omegaubd(discontubd) = NaN;
[~,discontufd] = max(abs(diff(omegaufd)));
omegaufd(discontufd) = NaN;
[~,discontdbd] = max(abs(diff(omegadbd)));
omegadbd(discontdbd) = NaN;
[~,discontdfd] = max(abs(diff(omegadfd)));
omegadfd(discontdfd) = NaN;
plot(endubdthetas,omegaubd,'b.')
plot(endubdthetas,omegaufd,'b-')
%plot(enddbdthetas,omegadbd,'r.')
plot(enddbdthetas,omegadfd,'r-')
scatter(endubdthetas(highlighti),omegaubd(highlighti),'bo')
scatter(endubdthetas(highlighti),omegadbd(highlighti),'ro')
scatter(endubdthetas(highlighti),omegaufd(highlighti),'bo')
scatter(endubdthetas(highlighti),omegadfd(highlighti),'ro')
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
ylim([-pi,pi])
legend('$\omega^\pm_{pre}$ w.r.t $m_1$',...
'$\omega^+_{post}$ w.r.t $m_1$',...%'$\omega^-_{pre}$ w.r.t $m_1$',...
'$\omega^-_{post}$ w.r.t $m_1$',...
'Highlighted + trajectory','Highlighted - trajectory',...
'Interpreter','Latex')
fplot(@(x)0,'k-')
%set(gca,'FontSize',20);
title('(c)','Interpreter','Latex')

%We plot the example trajectory in the inertial frame, using as our zero
%time the close encounter time (which is equal to the negation of the
%backwards detection circle intercept time)
[~,highlightceiu] = min(vecnorm(highlightptsu(1:2,:) - [1 - mu, 0]'));
[~,highlightceid] = min(vecnorm(highlightptsd(1:2,:) - [1 - mu, 0]'));
highlightcetu = highlightstdtimesu(highlightceiu);
highlightcetd = highlightstdtimesu(highlightceid);
%We have to use a for loop because annoyingly r2i doesn't work with
%multiple times at once
for i = 1:numel(highlightstdtimesu)
    highlightptsiu(:,i) = r2i(highlightstdtimesu(i) - highlightcetu,...
                            highlightptsu(:,i));
end
for i = 1:numel(highlightstdtimesd)
    highlightptsid(:,i) = r2i(highlightstdtimesd(i) - highlightcetd,...
                            highlightptsd(:,i));
end
subplot(2,2,4);
cplot(highlightptsiu,c,'b-')
hold on;
cplot(highlightptsid,c,'r-')
cplot(highlightptsiu(:,1),c,'ko')
cplot(highlightptsid(:,1),c,'ko')
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
legend('Highlighted + trajectory',...
       'Highlighted - trajectory',...
       'Backwards initial condition','Interpreter','Latex')
axis square;
%set(gca,'FontSize',20);
set(gca,'YAxisLocation','right');
title('(d)','Interpreter','Latex')