function len = bezier3_length(p0,p1,p2,s)
%   give bezier3_length under condition:
%   
%   p0 p1 p1 p2 是控制点
%
%   p0-p1  和  p1-p2 的距离相同

p0=p0(:);
p1=p1(:);
p2=p2(:);

dis = norm(p1-p0);
theta = acos((p2-p1)'*(p1-p0)/dis/dis);

x2 = 1+cos(theta);
y2 = sin(theta);

p0 = [0,0,0]';
p1 = [1,0,0]';
p2 = [x2,y2,0]';

norm(p2)

ds=0.0001;
s=0:ds:1;
[p,dp] = bezier3(p0,p1,p1,p2,s);
len = sum(sqrt(sum(diff(p')'.*diff(p')')))*dis;
len = sum(sqrt(sum(dp.*dp)))*ds*dis

dp1 = dp;

dl = 3*sqrt(2*x2*s.^4 - 4*x2*s.^3 + (2*x2 + 4)*s.^2 - 4*s + 1);
% integral 3*sqrt(2*x2*s.^4 - 4*x2*s.^3 + (2*x2 + 4)*s.^2 - 4*s + 1)
%   0->1

% check length 5
y = s.^2-s;
dl = 3*sqrt(2)*sqrt(x2*y.^2+2*y+0.5);
len5 = sum(dl)*ds*dis

% check length 6
y=-0.24999:0.00001:0;
dl = 3*sqrt(2)*0.5*sqrt(x2*y.^2+2*y+0.5)./sqrt(y+0.25);
len6 = sum(dl)*0.00001*dis*2

% check length 7
z=0.00001:0.00001:0.25;
dl = 1.5*sqrt(2)*sqrt(x2)*sqrt(z + (2/x2 - 1/2) + 1/16./z);
len7 = sum(dl)*0.00001*dis*2

% check length 8
z=0.00001:0.00001:0.25;
dl = 1.5*sqrt(2)*sqrt(x2)*sqrt(z + 1/16./z + (2/x2 - 1/2));
len7 = sum(dl)*0.00001*dis*2

sum(dp.*dp);
x0 = p0(1);
y0 = p0(2);
z0 = p0(3);

x1 = p1(1);
y1 = p1(2);
z1 = p1(3);

x2 = p2(1);
y2 = p2(2);
z2 = p2(3);

xa = 2*(x0+x2-2*x1);
ya = 2*(y0+y2-2*y1);
za = 2*(z0+z2-2*z1);

xb = 2*(x1-x0);
yb = 2*(y1-y0);
zb = 2*(z1-z0);

% dx = xa*s+xb;
% dy = ya*s+yb;
% dz = za*s+zb;
% len = sqrt(dx.*dx + dy.*dy + dz.*dz);

A = xa*xa + ya*ya + za*za;
B = 2*(xa*xb + ya*yb + za*zb);
C = xb*xb + yb*yb + zb*zb;
% len = sqrt(A*s^2 + B*s + C);

D=B/A;
E=C/A;
% len = sqrt(A) * integral 0->1 sqrt(s^2 + D*s + E) ds

F=D/2;
G=E-D*D/4;
len = sqrt(A)*sqrt((s+F).*(s+F) + G);

% len = sqrt(A) * integral F->1+F sqrt(s^2 + G) ds

t=s+F;
sum(sqrt(A)*sqrt(t.*t + G))*0.0001;

% integral of sqrt(x^2+a^2)  is:
%
% 0.5 * x * sqrt(x^2+a^2) + 0.5 * a^2 * ln(x+sqrt(a^2+x^2))
%


v=F;
lower = sqrt(A)*(0.5*v*sqrt(G+v*v)+0.5*G*log(v+sqrt(G+v*v)));
v=1+F;
upper = sqrt(A)*(0.5*v*sqrt(G+v*v)+0.5*G*log(v+sqrt(G+v*v)));

len = (upper - lower);

end








