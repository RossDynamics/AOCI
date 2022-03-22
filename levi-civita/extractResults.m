%Plots the forward and backwards interceptions for the ith matched set of
%+ and - initial conditions, where the initial conditions have been sorted
%according to the theta value of the up initial conditions's backward
%detection intercepts. Run fixedRadiusAnalysis first.

i = 100;

figure; 
hold on;

%c = cs(c,'s.o.v.dmode','position');

%We also plot the circle of up forward intercepts, just as a reference
%circle.
cplot(endufdsorted,c,'.b')

cplot(endubdsorted(:,i),c,'ob')
cplot(enddbdsorted(:,matchindex(i)),c,'or')
cplot(endufdsorted(:,i),c,'ob')
cplot(enddfdsorted(:,matchindex(i)),c,'or')

timesunfd = linspace(solunfdEnd{endubdi(i)}.x(1),solunfdEnd{endubdi(i)}.x(end),50000);
    [~,trajunfdpts] = lc2standard(timesunfd,deval(timesunfd,solunfdEnd{endubdi(i)}),...
                                        mu);
    cplot(trajunfdpts,c,'b')

timesdnfd = linspace(soldnfdEnd{enddbdi(matchindex(i))}.x(1),soldnfdEnd{enddbdi(matchindex(i))}.x(end),50000);
    [~,trajdnfdpts] = lc2standard(timesdnfd,deval(timesdnfd,soldnfdEnd{enddbdi(matchindex(i))}),...
                                        mu);
    cplot(trajdnfdpts,c,'r')

timesunbd = linspace(solunbdEnd{endubdi(i)}.x(1),solunbdEnd{endubdi(i)}.x(end),1000);
    [~,trajunbdpts] = lc2standard(timesunbd,deval(timesunbd,solunbdEnd{endubdi(i)}),...
                                        mu);
    cplot(trajunbdpts,c,'b')

timesdnbd = linspace(soldnbdEnd{enddbdi(matchindex(i))}.x(1),soldnbdEnd{enddbdi(matchindex(i))}.x(end),1000);
    [~,trajdnbdpts] = lc2standard(timesdnbd,deval(timesdnbd,soldnbdEnd{enddbdi(matchindex(i))}),...
                                        mu);
    cplot(trajdnbdpts,c,'r')