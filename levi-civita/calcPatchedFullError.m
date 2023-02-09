%Calculates and plots the maximum error between the patched conics orbital
%elements and the + and - trajectories (see fixedRadiusAnalysis and
%calcFlybyEndpoints) for different parameter values in order to ensure
%convergence. Whenever this script is run, all initializations of the
%variables rpd and scaling *must* be commented out in order to ensure that
%this script is able to pass those values to fixedRadiusAnalysis (yes, I
%know it would make more sense to make fixedRadiusAnalysis a function, but
%if I did that I'd be trapping all of those useful variables inside
%function scope unless I moved them all into the global workspace, and that
%is a theoretically interesting idea but I think it'd be unnecessarily
%complex to implement in practice)

clear variables;

c = Sun_Jupiter_LCR_PCR3BP_Context;

%The Sun-Jupiter distance in km.
L = 7.784e8;
%The radius of Jupiter in km.
Jrkm = 71492;
%The radius of Jupiter in denormalized units.
Jr = Jrkm / L

%The range of values to iterate over
valrangescaling = linspace(8,1,8);
valrangerpd = linspace(20, 400, 8);
%This switch specifies the variable this function should iterate over
mode = 'both'

errorenergypreplus = [];
errorenergypreminus = [];
errorenergypostplus = [];
errorenergypostminus = [];
erroromegapreplus = [];
erroromegapreminus = [];
erroromegapostplus = [];
erroromegapostminus = [];
curves = {};

%We use iterationvariable instead of i intentionally. Because of the
%horrible scoping stuff we're doing, fixedRadiusAnalysis will overwrite i,
%and so we need a variable name that it's not using.
for iterationvariable = 1:numel(valrangescaling)
    disp('===========================================')
    disp(iterationvariable)
    disp('===========================================')
    switch mode
        case 'scaling'
            %The default value for rpd is Jr/180
            rpd = Jr/180;
            scaling = valrangescaling(iterationvariable);
        case 'rpd'
            %The default value for scaling is 1
            rpd = Jr / valrangerpd(iterationvariable);
            scaling = 0.5;
        case 'both'
            scaling = valrangescaling(iterationvariable);
            rpd = Jr / valrangerpd(iterationvariable);
    end
    %We get the curves
    fixedRadiusAnalysis
    curves = calcFlybyEndpoints(detectDistStd,rpd,...
                                minusthetadotlist(100),n,c)
    close all;
    
    %and then calculate the max errors
    errorenergypreplus(iterationvariable) = ...
        max(abs(curves.energiespre - energiesubd));
    errorenergypreminus(iterationvariable) = ...
        max(abs(curves.energiespre - energiesdbd));    
    errorenergypostplus(iterationvariable) = ...
        max(abs(curves.energiespost - energiesufd));
    errorenergypostminus(iterationvariable) = ...
        max(abs(curves.energiespost - energiesdfd));
    
    %Because omega is an angle, we have to account for distance measuring
    %on a circle
    erroromegapreplus(iterationvariable) = ...
        max(abs(wrapToPi(curves.omegapre - omegaubd)));
    erroromegapreminus(iterationvariable) = ...
        max(abs(wrapToPi(curves.omegapre - omegadbd)));    
    erroromegapostplus(iterationvariable) = ...
        max(abs(wrapToPi(curves.omegapost - omegaufd)));
    erroromegapostminus(iterationvariable) = ...
        max(abs(wrapToPi(curves.omegapost - omegadfd)));

end

figure
hold on
plot(valrangescaling,errorenergypreplus,'bo-')
plot(valrangescaling,errorenergypreminus,'ro-')
plot(valrangescaling,errorenergypostplus,'b*-')
plot(valrangescaling,errorenergypostminus,'r*-')

legend('$E^+_{pre}$ error','$E^-_{pre}$ error','$E^+_{post}$ error',...
       '$E^-_{post}$ error','Interpreter','Latex')
xlabel('$r_d/r_h$','Interpreter','Latex')
ylabel('units','Interpreter','Latex')
set(gca,'Yscale','log')

figure
hold on
plot(valrangescaling,erroromegapreplus,'bo-')
plot(valrangescaling,erroromegapreminus,'ro-')
plot(valrangescaling,erroromegapostplus,'b*-')
plot(valrangescaling,erroromegapostminus,'r*-')


legend('$\omega^+_{pre}$ error','$\omega^-_{pre}$ error',...
       '$\omega^+_{post}$ error','$\omega^-_{post}$ error',...
       'Interpreter','Latex')
xlabel('$\frac{r_d}{r_h}$','Interpreter','Latex')
ylabel('rad','Interpreter','Latex')
set(gca,'Yscale','log')