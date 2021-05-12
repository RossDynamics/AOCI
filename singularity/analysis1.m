%Be sure to have your AOCI, CMDS, and FILM folders on your path before
%running this script. Also, be sure to have an active ER3BP context object.

% hold on
% l = pi/3;
% a = 10.9004/5.2;
% e = 0.8;
% findinit
% integplot(linspace(0,9.1765*2*pi,100000),[x0 0]',c)
% axis equal
% disp(1)
% 
% a = (10.9004+1e-4)/5.2;
% findinit
% integplot(linspace(0,11*2*pi,100000),[x0 0]',c)
% axis equal
% 
% disp(2)
% 
% a = (10.9004-1e-4)/5.2;
% findinit
% integplot(linspace(0,11*2*pi,100000),[x0 0]',c)
% axis equal

%%
hold on
l = pi/2;
e = 0.9;
a = 12.9457/5.2;
findinit
integplot(linspace(0,8.061*2*pi,100000),[x0 0]',c)
axis equal
disp(1)

a = 12.9456/5.2;
findinit
integplot(linspace(0,11*2*pi,100000),[x0 0]',c)
axis equal
disp(2)

a = 12.9458/5.2;
findinit
integplot(linspace(0,11*2*pi,100000),[x0 0]',c)
axis equal
disp(3)