%Analyzes where trajectories near the arches of chaos go.

%We nondimensionalize a.
a = 13 / 5.2044;

numConds = 5;

%lowere = 0.4;
%uppere = 0.9;

lowere = 0.58;
uppere = 0.63;

%lowere = 0.6158;
%uppere = 0.616;

c = getDichotomyIcsCMDS(numConds,lowere,uppere,a,c);

ics = [cg(c,'aoci.di.ics')
       zeros(1,size(cg(c,'aoci.di.ics'),2))]
   
   
c = startCaching(c);

hold on

%We plot the location of L4
integplot(linspace(0,4*pi,1000),[0.5 - cg(c,'p.mu') sqrt(3)/2 0 0 0].',...
                   c,'ob');

%We plot the line between m_1 and L4
fplot(@(x)sqrt(3)*x+sqrt(3)*cg(c,'p.mu'),'b')

% integanimate(linspace(0,4*pi,100),ics,c)

for i = 1:size(ics,2)
    integplot(linspace(0,6*pi,1000),ics(:,i),c,'k');
    %integplot(linspace(0,-3*pi,1000),ics(:,i),c,'k');
    disp(i)
end

c = stopCaching(c);