%input (a,e)
%output (x,y,xd,yd)

G=sqrt(a*(1-e^2)) ;
L=sqrt(a);

T=1/a+2*G ; % Tisserand
E=-0.5*T  ; % 3-body energy

% % BELOW assumes th=pi/3 and g_bar = 0
% th=pi/3;
% r=G^2/(1+e/2) ;
% 
% rd=-sqrt( -G^2/r^2 + 2/r - 1/L^2 ) ;
% thd=G/r^2 - 1 ;
% 
% l = acos( (1-r/a)/e ) - r*rd/L ;

% BELOW assumes mean anomoly, l=pi/3, and g_bar = 0
f = @(r) -G^2/r^2 + 2/r - 1/L^2 - a/r^2*(  acos( (1-r/a)/e ) - l)^2 ;
r=fzero(@(r) f(r),1) ;
rd=-sqrt( -G^2/r^2 + 2/r - 1/L^2 ) ;

thd=G/r^2 - 1 ;

% l = acos( (1-r/a)/e ) - r*rd/L ; % EQUALS pi/3 when negative rd taken

th = acos( (G^2/r - 1)/e ) ;

x=r*cos(th); y=r*sin(th);

dum = [x/r,-y; y/r,x]*[rd;thd] ;

xd=dum(1); yd=dum(2);

x0=[x y xd yd];






