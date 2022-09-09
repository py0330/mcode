function len = bezier2_length(p0,p1,p2,s)
%UNTITLED4 此处提供此函数的摘要
%   此处提供详细说明
p0=p0(:);
p1=p1(:);
p2=p2(:);


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
% len = sqrt(A) * integral F->1+F sqrt(s^2 + G) ds


% integral of sqrt(x^2+a^2)  is:
%
% 0.5 * x * sqrt(x^2+a^2) + 0.5 * a^2 * ln(x+sqrt(a^2+x^2))
%
v=F;
lower = sqrt(A)*(0.5*v*sqrt(G+v*v)+0.5*G*log(v+sqrt(G+v*v)));
v=s+F;
upper = sqrt(A)*(0.5*v*sqrt(G+v*v)+0.5*G*log(v+sqrt(G+v*v)));

len = (upper - lower);

end








