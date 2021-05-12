function Xin = rot2iner(X,T,M,del)
% Xin=rot2iner(X,T,M,del);
%
% Transformation from M body-centered rotating (nondim.) coordinates to
% inertial (nondim.) coordinates
%
% M = 2 : smaller mass(M2) centered inertial coordinates
% M = 1 : LARGER  mass(M1) centered inertial coordinates
% M = 0 : center-of-mass   centered inertial coordinates
%
% X   = state 4- or 6-vector matrix (in nondim. rotating coordinates)
% T   = nondimensional CR3BP time
% del = the rotation offset at time t=0 (in radians)
%
% Xrot = inv(B)*(Xin - A)
% see p.6 of Cassall(1996)
%
% NOTE: nondim. means nondimensional CR3BP units where 
%	sum of primaries' mass = 1
%	constant distance between primaries = 1
%	period of primaries' orbit = 2*pi
% 
%-----------------------------------------------------------------------
% CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
% with the LARGER MASS, M1 to the left of the origin at (-mu,0)
% and the smaller mass, M2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%       (rotating coords)
%
%                 L4
% -L3------M1--+-----L1--M2--L2-
%                 L5
%
% Shane Ross (revised 7.15.97)
%
% global mu
global mu

if nargin<=3, del=0; if nargin<=2, M=1; end; end; % Defaults

if     M==1, d= -mu; 	% distance to LARGER  primary in CR3BP (M1)
elseif M==2, d=1-mu;	% distance to smaller primary in CR3BP (M2)
elseif M==0, d=0;	% center-of-mass is the origin
end

[m,N]=size(X); 
if N==2, X=[X(:,1:2) zeros(m,4)]; X=expand(X); end
if N==4, XX=zeros(m,6); XX(:,1:2)=X(:,1:2); XX(:,4:5)=X(:,3:4); X=XX; end

c=cos(T + del); s=sin(T + del);

Xin(1:m,1)= c.*(X(:,1)-d) - s.*X(:,2);
Xin(1:m,2)= s.*(X(:,1)-d) + c.*X(:,2);
Xin(1:m,3)=X(:,3);
Xin(1:m,4)=-s.*(X(:,1)+X(:,5)-d) + c.*(X(:,4) - X(:,2));
Xin(1:m,5)= c.*(X(:,1)+X(:,5)-d) + s.*(X(:,4) - X(:,2));
Xin(1:m,6)=X(:,6);

if N==4, Xin=[Xin(:,1:2) Xin(:,4:5)]; end 
