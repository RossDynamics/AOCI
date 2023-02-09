function curves = calcFlybyEndpoints(rd,rce,thetadot,thetaDensity,c)
%Uses an analytical approach to plot the m1 orbital elements flyby of m_2
%(with close encounter radius rce) at detection radius rd as a function of
%close approach theta, where thetaDensity is the number of theta values to
%generate between 0 and 2*pi. At close encounter (r = rce) the spacecraft
%must have theta-direction velocity thetadot. c should be the current
%context object. Run fixedRadiusAnalysis first and be sure that the plot
%showing the change in Keplerian energy and argument of perigee is active.
%Returns a struct containing some of the relevant computed patched conics
%curves in display ordering.

%The two-body mu for m2 is the same as the three-body mu 
mu = cg(c,'p.mu');

v = rce * thetadot;
energy = 1/2*v^2 - mu/rce;

a = -mu / (2*energy);
e = rce*v^2/mu - 1;

p = a*(1-e^2);

%We calculate the true anomalies at each endpoint. nudpost is the
%post-close encounter true anomaly, and nudpre is the pre-close encounter
%true anomaly.

nupre = acos((p - rd)/(e*rd));
nupost = 2*pi - nupre;

%We now get the position and velocity vectors
rvecpost = [rd*cos(nupost)
             rd*sin(nupost)];
rvecpre = [rd*cos(nupre)
             rd*sin(nupre)];
         
vvecpost = sqrt(mu/p)*[-sin(nupost)
                        e+cos(nupost)];
vvecpre = sqrt(mu/p)*[-sin(nupre)
                        e+cos(nupre)];
                    
%We now attempt to determine the times at which the detection radius is
%encountered. We set the time of perigee passage equal to 0.

Epost = 2*atan(sqrt((1-e)/(1+e))*tan(nupost/2));
Epre = 2*atan(sqrt((1-e)/(1+e))*tan(nupre/2));

tpost = (Epost-e*sin(Epost))/sqrt(mu/a^3);
tpre = (Epre-e*sin(Epre))/sqrt(mu/a^3);

%We build out arrays of the state vectors in order to prepare to rotate
%them                    
statepostarray = [rvecpost
                  vvecpost]*ones(1,thetaDensity);
                 
stateprearray = [rvecpre
                 vvecpre]*ones(1,thetaDensity);

%We define a list of thetas and the rotation matrix that we need
thetas = linspace(0,2*pi,thetaDensity);

rotationth = @(theta)[cos(theta) -sin(theta) 0           0
                    sin(theta)  cos(theta) 0           0
                    0           0          cos(theta) -sin(theta)
                    0           0          sin(theta)  cos(theta)];             

%We now rotate according to each theta and then displace by the position
%and velocity of m2
for i = 1:thetaDensity
    statepostarrayrot(:,i) = rotationth(thetas(i)+tpost)*statepostarray(:,i);
    stateprearrayrot(:,i) = rotationth(thetas(i)+tpre)*stateprearray(:,i);
                         
    statepostarrayrf(:,i) = statepostarrayrot(:,i) + [cos(tpost)
                                                      sin(tpost)
                                                        -sin(tpost)
                                                       cos(tpost)];
                                                   
    stateprearrayrf(:,i) =  stateprearrayrot(:,i) + [cos(tpre)
                                                      sin(tpre)
                                                        -sin(tpre)
                                                       cos(tpre)];
end

%Now, we switch from velocity coordinates to momentum coordinates
% statepostarraymom = [statepostarrayrot(1,:)
%                      statepostarrayrot(2,:)
%                      statepostarrayrot(3,:)-statepostarrayrot(2,:)
%                      statepostarrayrot(4,:)+statepostarrayrot(1,:)];
% 
% stateprearraymom = [stateprearrayrot(1,:)
%                      stateprearrayrot(2,:)
%                      stateprearrayrot(3,:)-stateprearrayrot(2,:)
%                      stateprearrayrot(4,:)+stateprearrayrot(1,:)];                 
% 
% [apost,epost,omegapost,energiespost] = getae(statepostarraymom,1,1,c);
% [apre,epre,omegapre,energiespre] = getae(stateprearraymom,1,1,c);

%We now convert to the inertial frame and then convert the origin to that
%of m1
% statepostarrayrf = rotationth(-tpost) * statepostarrayrot
% 
% statepostarrayrf = r2i(tpost,statepostarrayrot) - [-mu*cos(tpost)
%                                                    -mu*sin(tpost)
%                                                    0
%                                                    0] - [0
%                                                          0
%                                                         (1-mu)*sin(tpost)
%                                                        -(1-mu)*cos(tpost)];
% stateprearrayrf = r2i(tpre,stateprearrayrot) - [-mu*cos(tpre)
%                                                 -mu*sin(tpre)
%                                                 0
%                                                 0] - [0
%                                                       0
%                                                       (1-mu)*sin(tpre)
%                                                      -(1-mu)*cos(tpre)];                                     

rm1pre = [stateprearrayrf(1:2,:); zeros(1,thetaDensity)];
vm1pre = [stateprearrayrf(3:4,:); zeros(1,thetaDensity)];
rm1post = [statepostarrayrf(1:2,:); zeros(1,thetaDensity)];
vm1post = [statepostarrayrf(3:4,:); zeros(1,thetaDensity)];

energiespost = 1/2*vecnorm(vm1post).^2 - (1-mu)*(vecnorm(rm1post)).^-1;
energiespre = 1/2*vecnorm(vm1pre).^2 - (1-mu)*(vecnorm(rm1pre)).^-1;

Hpre = cross(rm1pre,vm1pre);
Hpost = cross(rm1post,vm1post);

evecpre = (cross(vm1pre,Hpre)-(1-mu).*rm1pre./vecnorm(rm1pre))./(1-mu);
evecpost= (cross(vm1post,Hpost)-(1-mu).*rm1post./vecnorm(rm1post))./(1-mu);

epre = vecnorm(evecpre);
epost = vecnorm(evecpost);

omegapre = (Hpre(3,:) < 0)*2*pi + ...
           (-1).^(Hpre(3,:) < 0) .* atan2(evecpre(2,:),evecpre(1,:));
omegapost =(Hpost(3,:) < 0)*2*pi + ...
           (-1).^(Hpost(3,:) < 0) .* atan2(evecpost(2,:),evecpost(1,:));

%We sort the pre-encounter detection radius angles and corresponding points
[thetapre,sortind] = sort(wrapToPi(thetas-nupre))
[thetapost,sortind] = sort(wrapToPi(thetas-nupost))

%We copy the existing figure from fixedRadiusAnalysis so that we can edit
%it without destroying the old one
copyobj(gcf,0)

%We plot comparisons of the Keplerian energy before the encounter and
%after the encounter
%figure;
subplot(2,2,1);
hold on;
plot(thetapre,energiespre(sortind),'g.')
plot(thetapre,energiespost(sortind),'g-')
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('$E^\pm_{pre}$ w.r.t $m_1$',...
       '$E^+_{post}$ w.r.t $m_1$',...
       '$E^-_{post}$ w.r.t $m_1$','Highlighted + trajectory',...
       'Highlighted - trajectory','Highlighted + trajectory',...
       'Highlighted - trajectory','0',...
      '$E_{pre}$ (patched conics)',...
      '$E_{post}$ (patched conics)',...
       'Interpreter','Latex')

%We remove all information about the sample trajectories
axischildren = allchild(gca);
for i = 1:numel(axischildren)
    if isa(axischildren(i), 'matlab.graphics.chart.primitive.Scatter')
        delete(axischildren(i))
    end
end
   
subplot(2,2,3);
hold on;
[~,discontpre] = max(abs(diff(omegapre)));
omegapre(discontpre) = NaN;
[~,discontpost] = max(abs(diff(omegapost)));
omegapost(discontpost) = NaN;
plot(thetapre,omegapre(sortind),'g.')
plot(thetapre,omegapost(sortind),'g-')
xlabel('$\theta^+_{pre}$','Interpreter','Latex')
xticks([-pi,-pi/2,0,pi/2,pi])
xticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'})
xlim([-pi,pi])
legend('$\omega^\pm_{pre}$ w.r.t $m_1$',...
'$\omega^+_{post}$ w.r.t $m_1$',...
'$\omega^-_{post}$ w.r.t $m_1$','Highlighted + trajectory',...
       'Highlighted - trajectory','Highlighted + trajectory',...
       'Highlighted - trajectory','0',...
'$\omega_{pre}$ (patched conics)',...
 '$\omega_{post}$ (patched conics)',...
       'Interpreter','Latex')
title('(b)','Interpreter','Latex')
   
%We remove all information about the sample trajectories
axischildren = allchild(gca);
for i = 1:numel(axischildren)
    if isa(axischildren(i), 'matlab.graphics.chart.primitive.Scatter')
        delete(axischildren(i))
    end
end   
   
%We remove the subplots with the sample trajectories
delete(subplot(2,2,2))
delete(subplot(2,2,4))

%We now reshape the remaining subplots in the figure to get rid of the
%"slots" from which subplots had been removed
ax = subplot(2,2,3);
subplot(2,1,2,ax);
ax = subplot(2,2,1);
subplot(2,1,1,ax);

%If the figure window is maximized, we unmaximize it so that we can resize
fig = gcf;
fig.WindowState = 'normal';

%For consistency, we resize the figure window:
set(gcf,'Position',[100,100,700,500])

%Hacky and awful coding practice; please forgive me for my sins
sgtitle([strrep(sprintf(...
        '$r_d \\approx %0.4f $ units, $r_{ce} \\approx %0.4e',rd,rce),...
        'e-0','\times 10^{-') '}$ units'],'Interpreter','Latex')

curves = struct;
curves.energiespre = energiespre(sortind);
curves.energiespost = energiespost(sortind);
curves.omegapre = omegapre(sortind);
curves.omegapost = omegapost(sortind);
            
end