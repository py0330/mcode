%% test bezier2
p0 = [0.1,0.2,0.3];
p1 = [0.8,-0.2,0.4];
p2 = [0.4,0.5,0.3];

ds = 0.0001;
s = 0:ds:1;
[p, dp] = bezier2(p0',p1',p2',s);

% 数值积分1：基于贝塞尔曲线计算
len1 = sum(sqrt(sum(diff(p')'.*diff(p')')))

% 数值积分2：基于贝塞尔曲线导数计算
len2 = sum(sqrt(sum(dp.*dp)))*ds

% 数值积分3：
len3 = bezier2_length(p0,p1,p2,1)
%% test bezier3
p0 = [0.1,0.2,0.3];
p1 = [0.8,-0.2,0.4];
p2 = [0.4,0.5,0.3];

ds = 0.0001;
s = 0:ds:1;
[p, dp] = bezier3(p0',p1',p1',p2',s);

% 数值积分1：基于贝塞尔曲线计算
len1 = sum(sqrt(sum(diff(p')'.*diff(p')')))

% 数值积分2：基于贝塞尔曲线导数计算
len2 = sum(sqrt(sum(dp.*dp)))*ds

% 数值积分3：
len3 = bezier3_length(p0,p1,p2,s)
%% test bezier3 arcs
p0 = [0.1,0.2,0.3];
p1 = [0.8,-0.2,0.4];
p2 = [0.4,0.5,0.3];

% theta = pi/2;
% dis=3.15;
% p0 = [0,0,0];
% p1 = [1,0,0]*3.15;
% p2 = [1 + cos(theta),sin(theta),0]*3.15;

s=0:0.001:1;
[p, dp] = bezier3(p0',p1',p1',p2',s);

[arc,darc] = bezier3_s2arc(p0,p1,p2,s);

hold on

% real darc by bezier3
plot(s,sqrt(sum(dp.*dp)))

% arc
plot(s,arc);

% darc / ds 
plot(s,darc);

% darc / ds : using diff
plot(s(2:end)-0.0005,diff(arc)/0.001);

% d^2 arc / ds^2 
plot(s(3:end)-0.0015,diff(diff(arc))/0.001/0.001);

%% test bezier3_arc2s
theta = pi/2;
dis=3.15;

p0 = [0.1,0.2,0.3];
p1 = [1,0,0.8]*3.15;
p2 = [1 + cos(theta),sin(theta),0]*3.15;
clc
arc = bezier3_s2arc(p0,p1,p2,0.6128);
s = bezier3_arc2s(p0,p1,p2,arc)

%%
p0 = [0.1,0.2,0.3];
p1 = [0.8,-0.2,0.4];
p2 = [0.4,0.5,0.3];
p3 = [0.5,0.6,0.7];

[p, dp] = bezier3(p0',p1',p2',p3',s);

%% deduce bezier3 
axis equal;
hold on
scatter3(p0(1), p0(2), p0(3))
scatter3(p1(1), p1(2), p1(3))
scatter3(p2(1), p2(2), p2(3))

plot3(p(1,:)',p(2,:)',p(3,:)');
plot3(p3(1,:)',p3(2,:)',p3(3,:)');


%% deduce bezier3 
clear
syms p0 p1 p2 p3 s

p2 = p1

b0 = p0*(1-s) + p1*s
b1 = p1*(1-s) + p2*s
b2 = p2*(1-s) + p3*s

c0 = b0*(1-s) + b1*s
c1 = b1*(1-s) + b2*s

d0 = collect(c0*(1-s) + c1*s,s)

dd0 = diff(d0,s)

% ex = expand(dd0*dd0)
%% deduce bezier3 
syms x2 s

x1=1;
y2 = sqrt(x1*x1 - (x2-x1)*(x2-x1))
p0 = [0;0;0]
p1 = [x1;0;0]
p2 = [x2;y2;0]


b0 = p0*(1-s) + p1*s
b1 = p1*(1-s) + p1*s
b2 = p1*(1-s) + p2*s

c0 = b0*(1-s) + b1*s
c1 = b1*(1-s) + b2*s

d0 = c0*(1-s) + c1*s


dd0 = diff(d0,s)

ex = expand(dd0.*dd0)

ex = ex(2)+ex(1)

% ex = expand(ex/9/x1)



collect(ex,s)
% factor(ex,s)

% eqn = 2*x2*s^3 - 4*x2*s^2 + (2*x2 + 4)*s - 4
% factor(eqn,s)
%%
syms s x2

y = s^2-s;
eqn = expand(x2*y^2+2*y+0.5)
collect(eqn,s)
%%
syms z x2

y = z-0.25;
eqn = expand(x2*y^2+2*y+0.5)
collect(eqn,z)
%%
syms y
s=0.5+sqrt(y+0.25)

expand(s*s-s)
%2*x2*s^4 - 4*x2*s^3 + (2*x2 + 4)*s^2 - 4*s + 1

%%
clear
syms theta s
% s=0.5
% theta = pi

p0 = [0,0,0];
p1 = [1,0,0];
p2 = [1,0,0];
p3 = [1 + cos(theta),sin(theta),0];


eqn = 3*p1 - 3*p0 - 3*(p0 - 3*p1 + 3*p2 - p3)*s.^2 + 2*(3*p0 - 6*p1 + 3*p2)*s;
eq1 = (3*cos(theta) + 3)*s^2 - 6*s + 3
eq2 = 3*s^2*sin(theta)

eq = eq1*eq1 + eq2*eq2
collect(eq,s)
sqrt(eq)


subs(sqrt(eq),theta,pi)
%%
clear
syms A B C D s

eqn = sqrt(A*s^4 + B*s^3 + C*s^2 + D*s + 1)
dif_dqn = diff(eqn,s)


subs(dif_dqn,s,0)
%%
clear;
s=0:0.00001:1;

hold on
axis equal
plot(s,(sin(s*pi) + 0.5*cos(s*2*pi) - 0.5)/pi)
%% 推导bezier曲率
clear
syms s;
syms x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3
p0 = [x0;y0;z0];
p1 = [x1;y1;z1];
p2 = [x2;y2;z2];
p3 = [x3;y3;z3];




% syms p0 p1 p2 p3 s
p = (3*p1 - p0 - 3*p2 + p3)*s.^3 + (3*p0 - 6*p1 + 3*p2)*s.^2 + (3*p1 - 3*p0)*s + p0;
dp = 3*p1 - 3*p0 - 3*(p0 - 3*p1 + 3*p2 - p3)*s.^2 + 2*(3*p0 - 6*p1 + 3*p2)*s;
ddp = 6*p0 - 12*p1 + 6*p2 - 2*s*(3*p0 - 9*p1 + 9*p2 - 3*p3);

darc = sqrt(dp(1)^2 + dp(2)^2 + dp(3)^2);
ds_over_darc = 1/darc;
ds2_over_darc2 = diff(ds_over_darc,s).*ds_over_darc

dp_over_darc = dp * ds_over_darc;
dp2_over_darc2 = ddp*ds_over_darc^2 + dp * ds2_over_darc2;

%%
s = 0
% ds_over_darc
1/(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2)^(1/2)

(2*(6*x0 - 12*x1 + 6*x2 - 2*s*(3*x0 - 9*x1 + 9*x2 - 3*x3))*((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1) + 2*(6*y0 - 12*y1 + 6*y2 - 2*s*(3*y0 - 9*y1 + 9*y2 - 3*y3))*((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1) + 2*(6*z0 - 12*z1 + 6*z2 - 2*s*(3*z0 - 9*z1 + 9*z2 - 3*z3))*((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1))/(2*(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2)^2)

%%

% x1 = x2;
% y1 = y2;
% z1 = z2;
s=0

syms xa ya za xb yb zb

x1 = x0 + xa;
y1 = y0 + ya;
z1 = z0 + za;

x2 = x1 + xb;
y2 = y1 + yb;
z2 = z1 + zb;

n1 = xa^2 + ya^2 + za^2;

dp2_over_darc2 = ...
[
(6*x0 - 12*x1 + 6*x2 - 2*s*(3*x0 - 9*x1 + 9*x2 - 3*x3))/(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2) - ((2*(6*x0 - 12*x1 + 6*x2 - 2*s*(3*x0 - 9*x1 + 9*x2 - 3*x3))*((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1) + 2*(6*y0 - 12*y1 + 6*y2 - 2*s*(3*y0 - 9*y1 + 9*y2 - 3*y3))*((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1) + 2*(6*z0 - 12*z1 + 6*z2 - 2*s*(3*z0 - 9*z1 + 9*z2 - 3*z3))*((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1))*((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1))/(2*(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2)^2)
(6*y0 - 12*y1 + 6*y2 - 2*s*(3*y0 - 9*y1 + 9*y2 - 3*y3))/(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2) - ((2*(6*x0 - 12*x1 + 6*x2 - 2*s*(3*x0 - 9*x1 + 9*x2 - 3*x3))*((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1) + 2*(6*y0 - 12*y1 + 6*y2 - 2*s*(3*y0 - 9*y1 + 9*y2 - 3*y3))*((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1) + 2*(6*z0 - 12*z1 + 6*z2 - 2*s*(3*z0 - 9*z1 + 9*z2 - 3*z3))*((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1))*((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1))/(2*(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2)^2)
(6*z0 - 12*z1 + 6*z2 - 2*s*(3*z0 - 9*z1 + 9*z2 - 3*z3))/(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2) - ((2*(6*x0 - 12*x1 + 6*x2 - 2*s*(3*x0 - 9*x1 + 9*x2 - 3*x3))*((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1) + 2*(6*y0 - 12*y1 + 6*y2 - 2*s*(3*y0 - 9*y1 + 9*y2 - 3*y3))*((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1) + 2*(6*z0 - 12*z1 + 6*z2 - 2*s*(3*z0 - 9*z1 + 9*z2 - 3*z3))*((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1))*((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1))/(2*(((3*x0 - 9*x1 + 9*x2 - 3*x3)*s^2 + (12*x1 - 6*x0 - 6*x2)*s + 3*x0 - 3*x1)^2 + ((3*y0 - 9*y1 + 9*y2 - 3*y3)*s^2 + (12*y1 - 6*y0 - 6*y2)*s + 3*y0 - 3*y1)^2 + ((3*z0 - 9*z1 + 9*z2 - 3*z3)*s^2 + (12*z1 - 6*z0 - 6*z2)*s + 3*z0 - 3*z1)^2)^2)
]

% syms n
% (6*xa)/(9*n^2) - (3*xa*(36*n^2))/(2*(9*n^2)^2)
dp2_over_darc2 - [
(3*xa*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*xa - 6*xb)/(9*n1)
(3*ya*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*ya - 6*yb)/(9*n1)
(3*za*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*za - 6*zb)/(9*n1)
]

syms n1

test = [
(3*xa*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*xa - 6*xb)/(9*n1)
(3*ya*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*ya - 6*yb)/(9*n1)
(3*za*(3*xa*(12*xa - 12*xb) + 3*ya*(12*ya - 12*yb) + 3*za*(12*za - 12*zb)))/(2*(9*n1)^2) - (6*za - 6*zb)/(9*n1)
]
expand(test)

((xb*ya^2) + (xb*za^2) - (xa*ya*yb) - (xa*za*zb))     /   (3*n1^2)    *  2



%%
subplot(1,2,1)

p0 = [2,0,0];
p1 = [2,2,0];
p2 = [1.0,2,0];
p3 = [0,2,0];

s = 0:0.001:1;
[p, dp] = bezier3(p0',p1',p2',p3',s);

axis equal;
hold on
scatter3(p0(1), p0(2), p0(3))
scatter3(p1(1), p1(2), p1(3))
scatter3(p2(1), p2(2), p2(3))
scatter3(p3(1), p3(2), p3(3))

plot3(p(1,:)',p(2,:)',p(3,:)');
% plot3(p3(1,:)',p3(2,:)',p3(3,:)');

p = 2*[cos(s*pi)',sin(s*pi)',zeros(length(s),1)]'
plot3(p(1,:)',p(2,:)',p(3,:)');
% 
% dp2_over_darc2

%%
subplot(1,2,1)
ratio = 1.5;

p0 = [2,0,0];
p1 = [2,1*ratio,0];
p2 = [2-0.75*ratio*ratio,2*ratio,0];
p3 = [4,2,2.5];

ds = 0.001;
s = 0:ds:1;
[p, dp] = bezier3(p0',p1',p2',p3',s);

dp_over_ds = diff(p')'/ds;
d2p_over_ds2 = diff(dp_over_ds')'/ds;

darc_over_ds = sqrt(sum(dp_over_ds.^2));
d2arc_over_ds2 = diff(darc_over_ds) / ds;

ds_over_darc = 1./darc_over_ds;
d2s_over_darc2 = -d2arc_over_ds2./darc_over_ds(1:end-1).^2;

dp_over_darc = dp_over_ds.*ds_over_darc;
d2p_over_darc2 = d2p_over_ds2.*ds_over_darc(1:end-1).^2 + dp_over_ds(:,1:end-1) .* d2s_over_darc2;
% ddp_over_darc = diff(dp_over_darc')'/ds./(0.5*darc_over_ds(1:end-2) + 0.5*darc_over_ds(2:end-1));

axis equal;
hold on
scatter3(p0(1), p0(2), p0(3))
scatter3(p1(1), p1(2), p1(3))
scatter3(p2(1), p2(2), p2(3))
scatter3(p3(1), p3(2), p3(3))

plot3(p(1,:)',p(2,:)',p(3,:)');

r = 2;
p2 = r*[cos(s*pi)',sin(s*pi)',zeros(length(s),1)]';
plot3(p2(1,:)',p2(2,:)',p2(3,:)');

dp2 = diff(p2')'/ds;
darc_over_ds2 = sqrt(sum(dp2.^2));
dp_over_darc2 = dp2./darc_over_ds2;
ddp_over_darc2 = diff(dp_over_darc2')'/ds./(0.5*darc_over_ds2(1:end-1) + 0.5*darc_over_ds2(2:end));

subplot(1,2,2)
% axis equal
hold on
% plot(s(2:end)-ds/2,dp2)
% plot(s(3:end)-ds*3/2,ddp_over_darc2)
plot(s(3:end)-ds*3/2,d2p_over_darc2)
% plot(s(2:end)-ds/2,ddp_over_darc/ds)
% plot(s(2:end)-ds/2,diff(p')'/ds)
% plot(s(3:end)-ds*3/2,diff(diff(p'))'/pi^2/r^2/ds/ds)
% plot(s(2:end)-ds/2, ddp_over_darc)

%%
[p, dp] = bezier3(p0',p1',p2',p3',s);
p_ = 2*[cos(s*pi)',sin(s*pi)',zeros(length(s),1)]';

subplot(1,2,1)
hold on
plot(diff(p(2,:)')/3)
plot(diff(p_(2,:)')/2/pi)

% subplot(1,2,2)
% hold on
% plot(diff(p(1,:)')/3)
% plot(diff(p_(1,:)')/2/pi)

