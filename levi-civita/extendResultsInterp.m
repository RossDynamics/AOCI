%Extends the solution to the ith down trajectory over a longer timespan. 
%Run fixedRadiusAnalysisInterp first.

i = mini-2;

newend = 2 * endtime;

%figure; 
hold on;

soldnbdEndext = odextend(solintdnfdEnd{i},[],newend);


timesdnfd = linspace(soldnbdEndext.x(1),soldnbdEndext.x(end),100000);
    [~,trajdnfdpts] = lc2standard(timesdnfd,deval(timesdnfd,soldnbdEndext),...
                                        mu);
    cplot(trajdnfdpts,c,'r');