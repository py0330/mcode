%% deduce s_blend_line_circle
syms p0 p1 s center theta
% line part
d0 = p0*(1-s)+p1*s;
d1 = p1;
d2 = p1;

e0 = d0*(1-s)+d1*s;
e1 = p1;

f = e0*(1-s) + e1 *s;


line_part = f;

expand(f)
expand(line_part)

collect(f,s)

collect(diff(diff(f,s)),s)

syms rx ry
diff(sin(s*s*s*theta).*ry + cos(s*s*s*theta).*rx, s)
diff(diff(sin(s*s*s*theta).*ry + cos(s*s*s*theta).*rx, s),s)

% 
% p2 = p1;
% 
% b0 = p0*(1-s) + p1*s
% b1 = p1*(1-s) + p2*s
% b2 = center + (sin((1-s)*theta).*(p2 - center) + sin(s*theta).*(p3 - center)) / sin(theta)
% 
% c0 = b0*(1-s) + b1*s
% c1 = b1*(1-s) + b2*s
% 
% d0 = collect(c0*(1-s) + c1*s,s)
% 
% dd0 = diff(d0,s)
% d2d0 = diff(dd0,s)

%% s_blend_line_line
clear
p0=[0.3;0.8;-0.3];
p1=[-0.5;0.3;-0.3];
p2=[0.3;-0.5;0.2];

ds = 0.01;
s = 0:ds:1;
[p,dp,ddp] = s_blend_line_line_bezier3(p0,p1,p2,s);

subplot(1,3,1)
hold on
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

% plot line
line = p0*s+p1*(1-s);
plot3(line(1,:),line(2,:),line(3,:))

% plot circle
line = p1*s+p2*(1-s);
plot3(line(1,:),line(2,:),line(3,:))

subplot(1,3,2)
hold on
plot(diff(p')/ds)
plot(dp')

subplot(1,3,3)
hold on
plot(diff(dp')/ds)
plot(ddp')

hold off
[ds_over_darc, d2s_over_darc2, darc, d2arc] = s_ds_over_darc(dp, ddp);

plot(darc)
hold on
plot(d2arc)

darc_compare = sqrt(sum(dp.^2));
plot(darc_compare,'.--')
plot(diff(darc_compare)/ds,'.--')


[arc_est, darc_est] = s_estimate_bezier3_arc(darc(1), d2arc(1), darc(end), d2arc(end), darc(0.5/ds + 1), s);

% A = [0,0,0,1;1,1,1,1;0,0,1,0;3,2,1,0];
% 
% coe = A^-1 * [darc(1), darc(end), d2arc(1), d2arc(end)]';
% a = coe(1)
% b = coe(2)
% c = coe(3)
% d = coe(4)
% darc_estimate = a.*s.^3 + b.*s.^2 + c.*s + d;


plot(darc_est,'.--')
plot((1.5:1:size(s,2))',diff(arc_est)/ds,'.-.')

% plot(ds_over_darc)
% hold on
% plot(d2s_over_darc2)
% hold on
% ds_compare = 1./sqrt(sum(dp.^2));
% % plot(ds_compare,'.')
% plot(diff(ds_compare)./(darc(1:end-1)'*ds),'.')



%% s_blend_line_circle
clear
theta = pi/18;
center = [0,0,0]';
radius = 5;
ax = [0,0,1]';

rkkk = cross([-1,0,0]',ax);
rkkk = rkkk/norm(rkkk);

p0 = [3.5,0.7,-2.8]';
p1 = center + radius * rkkk;

p0 = cross(p1-center,ax);
p0 = p0/norm(p0)*theta*radius;
p0 = p1 + p0;


ds = 0.01;
s = 0:ds:1;
[p,dp,ddp] = s_blend_line_circle_bezier3(p0,p1,center,ax,theta,s);

subplot(1,3,1)
hold on
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

scatter3(p1(1),p1(2),p1(3))
scatter3(p0(1),p0(2),p0(3))

% plot line
line = p0*s+p1*(1-s);
plot3(line(1,:),line(2,:),line(3,:))

% plot circle
rx = p1 - center;
ry = cross(ax, rx);
circle = center + cos(theta.*s).*rx+sin(theta.*s).*ry;
plot3(circle(1,:),circle(2,:),circle(3,:))

subplot(1,3,2)
hold on
plot(diff(p')/ds)
plot(dp')

subplot(1,3,3)
hold on
plot(diff(dp')/ds)
plot(ddp')

hold off
[ds_over_darc, d2s_over_darc2, darc, d2arc] = s_ds_over_darc(dp, ddp);

plot(darc)
hold on
plot(d2arc)

darc_compare = sqrt(sum(dp.^2));
plot(darc_compare,'.--')
plot(diff(darc_compare)/ds,'.--')


[arc_est, darc_est] = s_estimate_bezier3_arc(darc(1), d2arc(1), darc(end), d2arc(end), darc(0.5/ds + 1), s);

% A = [0,0,0,1;1,1,1,1;0,0,1,0;3,2,1,0];
% 
% coe = A^-1 * [darc(1), darc(end), d2arc(1), d2arc(end)]';
% a = coe(1)
% b = coe(2)
% c = coe(3)
% d = coe(4)
% darc_estimate = a.*s.^3 + b.*s.^2 + c.*s + d;


plot(darc_est,'.--')
plot((1.5:1:size(s,2))',diff(arc_est)/ds,'.-.')

% plot(ds_over_darc)
% hold on
% plot(d2s_over_darc2)
% hold on
% ds_compare = 1./sqrt(sum(dp.^2));
% % plot(ds_compare,'.')
% plot(diff(ds_compare)./(darc(1:end-1)'*ds),'.')

%% s_blend_circle_circle
clc
clear

p1 = [0.8;0.3;-1.5];

theta1 = pi/6;
c1 = [-1.5;0.8;0.2];
ax1 = cross(c1-p1,[0.5;0.3;0.8]);
ax1 = ax1/norm(ax1);

theta2 = pi*2/3;
c2 = [1.2;0.4;-0.3];
ax2 = cross(c2-p1,[0;1;0]);
ax2 = ax2/norm(ax2);

ds = 0.01;
s = 0:ds:1;
[p,dp,ddp] = s_blend_circle_circle_bezier3(p1,c1,ax1,theta1,c2,ax2,theta2,s);

subplot(1,3,1)
hold on
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

plot3(p1(1,:),p1(2,:),p1(3,:), '*')

rx = p1 - c1;
ry = -cross(ax1, rx);
circle = c1 + cos(theta1.*s).*rx+sin(theta1.*s).*ry;
plot3(circle(1,:),circle(2,:),circle(3,:))

rx = p1 - c2;
ry = cross(ax2, rx);
circle = c2 + cos(theta2.*s).*rx+sin(theta2.*s).*ry;
plot3(circle(1,:),circle(2,:),circle(3,:))


figure1 = plot3(p(1,1),p(2,1),p(3,1), '.','MarkerFaceColor','r','MarkerSize',10)

% for k = 1:size(p,2)
%     figure1.XData = p(1,k);
%     figure1.YData = p(2,k);
%     figure1.ZData = p(3,k);
%     drawnow
%     pause(0.1);
% end

subplot(1,3,2)
hold on
plot(diff(p')/ds)
plot(dp')

subplot(1,3,3)
hold on
plot(diff(dp')/ds)
plot(ddp')


%% s_blend_circle_circle
clc
clear

% center = [0.57345962137149420;0.22504954769300603;0.23562075307425620]
% p0 = [0.60198384418161333;0.20742573893289451;0.24832503880025239]
% axis = [0.0019395681385318642;-0.58269493919875459;-0.81268865250327915]
% radius = 0.035855666190186097

p1 = [0.8;0.3;-1.5];

theta1 = 0.01/0.035855666190186097;
c1 = [-1.5;0.8;0.2];
ax1 = cross(c1-p1,[0.5;0.3;0.8]);
ax1 = ax1/norm(ax1);

theta2 = pi*2/3;
c2 = [0.58403509793717567;0.22044018761897119;-0.3];
ax2 = cross(c2-p1,[0;1;0]);
ax2 = ax2/norm(ax2);

ds = 0.01;
s = 0:ds:1;
[p,dp,ddp] = s_blend_circle_circle_bezier3(p1,c1,ax1,theta1,c2,ax2,theta2,s);

subplot(1,3,1)
hold on
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

plot3(p1(1,:),p1(2,:),p1(3,:), '*')

rx = p1 - c1;
ry = -cross(ax1, rx);
circle = c1 + cos(theta1.*s).*rx+sin(theta1.*s).*ry;
plot3(circle(1,:),circle(2,:),circle(3,:))

rx = p1 - c2;
ry = cross(ax2, rx);
circle = c2 + cos(theta2.*s).*rx+sin(theta2.*s).*ry;
plot3(circle(1,:),circle(2,:),circle(3,:))


figure1 = plot3(p(1,1),p(2,1),p(3,1), '.','MarkerFaceColor','r','MarkerSize',10)

% for k = 1:size(p,2)
%     figure1.XData = p(1,k);
%     figure1.YData = p(2,k);
%     figure1.ZData = p(3,k);
%     drawnow
%     pause(0.1);
% end

subplot(1,3,2)
hold on
plot(diff(p')/ds)
plot(dp')

subplot(1,3,3)
hold on
plot(diff(dp')/ds)
plot(ddp')





%% s_blend_line_circle arc 




%% s_blend_line_line_reverse
clear
theta = 0.1;

p0=[1;0;0];
p1=[0;0;0];
p2=[cos(theta);sin(theta);0];

ds = 0.01;
s = 0:ds:1;
[p,dp,ddp] = s_blend_line_line_bezier3(p0,p1,p2,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 1
subplot(1,3,1)
hold on
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

% plot line 1
line = p0*s+p1*(1-s);
plot3(line(1,:),line(2,:),line(3,:))

% plot line 2
line = p1*s+p2*(1-s);
plot3(line(1,:),line(2,:),line(3,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 2
subplot(1,3,2)
hold on
plot(diff(p')/ds)
plot(dp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3
subplot(1,3,3)
hold on
plot(diff(dp')/ds)
plot(ddp')

hold off
[ds_over_darc, d2s_over_darc2, darc, d2arc] = s_ds_over_darc(dp, ddp);

plot(darc)
hold on
plot(d2arc)

darc_compare = sqrt(sum(dp.^2));
plot(darc_compare,'.--')
plot(diff(darc_compare)/ds,'.--')


[arc_est, darc_est, ddarc_est] = s_estimate_bezier3_arc(darc(1), d2arc(1), darc(end), d2arc(end), darc(0.5/ds + 1), s);

% A = [0,0,0,1;1,1,1,1;0,0,1,0;3,2,1,0];
% 
% coe = A^-1 * [darc(1), darc(end), d2arc(1), d2arc(end)]';
% a = coe(1)
% b = coe(2)
% c = coe(3)
% d = coe(4)
% darc_estimate = a.*s.^3 + b.*s.^2 + c.*s + d;


plot(darc_est,'.--')
plot(ddarc_est,'.--')
plot((1.5:1:size(s,2))',diff(arc_est)/ds,'.-.')

% c = theta*theta/8/pi;
% a = norm(ddp(:,end)) / atan(1.0/c/2);
% 
% 
% x = s-0.5;
% plot(atan(x/c) *a)

% plot(ds_over_darc)
% hold on
% plot(d2s_over_darc2)
% hold on
% ds_compare = 1./sqrt(sum(dp.^2));
% % plot(ds_compare,'.')
% plot(diff(ds_compare)./(darc(1:end-1)'*ds),'.')
%%
clear
m = 3;
n = 4;
ds = 0.001;
for i=1:m*n
%     i
%     m
%     n
    theta = linspace(0,pi,m*n);
    theta = theta(i);
    
    r = 2.0;
    p0=[r;0;0];
    p1=[0;0;0];
    p2=[r*cos(theta);r*sin(theta);0];

    s = 0:ds:1;
    [p,dp,ddp] = s_blend_line_line_bezier3(p0,p1,p2,s);
    subplot(m,n,i)
    hold on;
    [ds_over_darc, d2s_over_darc2, darc, d2arc] = s_ds_over_darc(dp, ddp);
    plot(s,darc);
    plot(s,d2arc);
    darc_compare = sqrt(sum(dp.^2));
    plot(s(1:end),darc_compare,'.--');
    plot(s(1:end-1)+ds/2,diff(darc_compare)/ds,'.--');
    

    theta = (darc(1/ds/2+1)/darc(end))*2*pi;
    c = theta^1.5*0.06;
    a = norm(ddp(:,end)) / atan(1.0/c/2);

    x = s-0.5;
    plot(s,atan(x/c) *a);

    [arc_est, darc_est, ddarc_est] = s_estimate_bezier3_arc(darc(1), d2arc(1), darc(end), d2arc(end), darc(0.5/ds + 1), s);
    plot(s(1:end-1)+ds/2,diff(arc_est)/ds,'--');
    plot(s,darc_est,'--');
    plot(s,ddarc_est,'--');
    title([num2str(theta)]);
    
%     darc_est(501)
%     darc(501) - darc_est(501)
    arc_est(1)
end
%%
% H = -(darc50 - ratio*(a/8 + b/4 + c/2 + d) + (j - h*i*log(h))*(ratio - 1))/(4*pi);
% X = (ratio - 1)*(j/2 - k + (h*i)/4)/(4*pi);

clear
syms a b c d h i j k lx at s ratio darc50
ddarcA = a.*s.^2.*3 + b.*s   .*2 + c                ;
darcA  = a.*s.^3    + b.*s.^2    + c.*s       + d   ;
arcA   = a.*s.^4./4 + b.*s.^3./3 + c.*s.^2./2 + d.*s;

x = s - 0.5;

ddarcB = at;
darcB  = x.*at - lx + j;
arcB   = k - at*h^2/2 + x.*x.*at/2 - x.*lx + x*i*h/2+ x*j;

darcA50 = 0.125.*a + 0.25.*b + 0.5.*c + d;
darcB50 = -log(h)*h*i + j;
darcE50 = darc50 - (ratio * darcA50   + (1-ratio)*darcB50);

arc   = ratio*arcA    + (1-ratio)*arcB   + s*darcE50/2 - sin(2*pi*s)/2*darcE50/2/pi;
darc  = ratio*darcA   + (1-ratio)*darcB  + (1-cos(2*pi*s))/2 * darcE50;
ddarc = ratio*ddarcA  + (1-ratio)*ddarcB + sin(2*pi*s)*pi*darcE50;

expand(arc)
D = darc50/2 - (ratio - 1)*(j + (h*i)/2) + d*ratio - (ratio*(a/8 + b/4 + c/2 + d))/2 + ((j - h*i*log(h))*(ratio - 1))/2;
collect(darc,s)
Y = d*ratio - (cos(2*pi*s)/2 - 1/2)*(darc50 - ratio*(a/8 + b/4 + c/2 + d) + (j - h*i*log(h))*(ratio - 1)) + (ratio - 1)*(at/2 - j + lx);
collect(Y,s)
D = darc50/2 - (ratio - 1)*(j + (h*i)/2) + d*ratio - (ratio*(a/8 + b/4 + c/2 + d))/2 + ((j - h*i*log(h))*(ratio - 1))/2;
J = -(darc50 - ratio*(a/8 + b/4 + c/2 + d) + (j - h*i*log(h))*(ratio - 1))/(4*pi);
d*ratio + - 1/2*J*4*pi + (ratio - 1)*( - j) - D


collect(ddarc,s)



%%
syms x h i k
f(x)=i*(...
        (x.*x+h*h)/2.*atan2(x, h)-h*x/2 ...
        - (x.*log(x.*x + h*h)+2*h*atan2(x, h)-2*x)*h/2 ...
        ) + k;
f(-0.5)

eqn = i*(...
        (x.*x+h*h)/2.*atan2(x, h)-h*x/2 ...
        - (x.*log(x.*x + h*h)+2*h*atan2(x, h)-2*x)*h/2 ...
        ) + k
expand(eqn)
%%
c = theta/8/pi;
a = 6 / atan(1.0/c/2);

x = s-0.5;

ddp = atan(x/c) *a;
dp = (x.*atan(x/c) - c/2*log(x.*x + c*c))*a;
p  = a*(...
        (x.*x+c*c)/2.*atan(x/c)-c*x/2 ...
        - (x.*log(x.*x + c*c)+2*c*atan(x/c)-2*x)*c/2 ...
        );


figure
hold on
plot(s,ddp)
plot(s(1:end-1)+ds/2, diff(dp)/ds)
plot(s(1:end-2)+ds, diff(diff(p))/ds/ds)





% dp = 

%% bezier 3 blend line circle
theta = pi/3;
center = [0;0;0];
radius = 1;

p0 = [3.5,0.7,-2.8]';
p1 = radius * [1,0,0]';
p2 = p1;
p3 = p2;
p4 = radius * [cos(theta),sin(theta),0]';

ds = 0.01;
s = 0:ds:1;

b0 = p0.*(1-s) + p1.*s;
b1 = p1.*(1-s) + p2.*s;
b2 = p2.*(1-s) + p3.*s;
b3 = center + (sin((1-s)*theta).*(p3 - center) + sin(s.*theta).*(p4 - center)) ./ sin(theta);

c0 = b0.*(1-s) + b1.*s;
c1 = b1.*(1-s) + b2.*s;
theta = acos(sum((b2 - center).*(b3 - center))./radius^2);
c2 = center + (sin((1-s).*theta).*(b2 - center) + sin(s.*theta).*(b3 - center)) ./ sin(theta);

d0 = c0.*(1-s) + c1.*s;
theta = acos(sum((c1 - center).*(c2 - center))./radius^2);
d1 = center + (sin((1-s).*theta).*(c1 - center) + sin(s.*theta).*(c2 - center)) ./ sin(theta);

e = d0 + d1 - p1;

hold on

% p = c2;
% plot3(p(1,:),p(2,:),p(3,:) ,'.')
% axis equal
% 
% p = b3;
% plot3(p(1,:),p(2,:),p(3,:) ,'.')
% axis equal

p = d1;
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

p = d0;
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

p = e;
plot3(p(1,:),p(2,:),p(3,:), '.--')
axis equal

d1(:,1)

%%
theta = pi/3;
center = [0;0;0];

p0 = [3.5,0.7,-2.8]';
p1 = [1,0,0]';
p2 = p1;
p3 = [cos(theta),sin(theta),0]';

ds = 0.0001;
s = 0:ds:1;

b0 = p0.*(1-s) + p1.*s
b1 = p1.*(1-s) + p2.*s
b2 = center + (sin((1-s)*theta).*(p2 - center) + sin(s.*theta).*(p3 - center)) ./ sin(theta)

c0 = b0.*(1-s) + b1.*s
c1 = center + (sin((1-s)*theta).*(b1 - center) + sin(s.*theta).*(b2 - center)) ./ sin(theta)

d0 = c0.*(1-s) + c1.*s
d0 = center + (sin((1-s)*theta).*(c0 - center) + sin(s.*theta).*(c1 - center)) ./ sin(theta)

dd0 = diff(d0')'
d2d0 = diff(dd0')'

dd0(:,1)/norm(dd0(:,1))
(p0-p1)/norm(p0-p1)
dd0(:,1)/norm(dd0(:,1)) + (p0-p1)/norm(p0-p1)


hold on

p = b2;
plot3(p(1,:),p(2,:),p(3,:))
axis equal

p = c1;
plot3(p(1,:),p(2,:),p(3,:))
axis equal

p = d0;
plot3(p(1,:),p(2,:),p(3,:))
axis equal

%%
s = 1
(p1 - p0)*s^3 + (center + 3*p0 - 4*p1 + (sin(theta*(s - 1))*(center - p1) - sin(s*theta)*(center - p3))/sin(theta))*s^2 + (3*p1 - 3*p0)*s + p0
%% diff bezier
s = 1;
center = 0;

syms x
x = 1
p1 = [x;0;0];
p3 = [x*cos(theta); x*sin(theta);0];

3*p1 - 3*p0 - 3*s^2*(p0 - p1) + 2*s*(center + 3*p0 - 4*p1 + (sin(theta*(s - 1))*(center - p1) - sin(s*theta)*(center - p3))/sin(theta)) - (s^2*(theta*cos(s*theta)*(center - p3) - theta*cos(theta*(s - 1))*(center - p1)))/sin(theta)
%% diff 2 bezier


2*center + 6*p0 - 8*p1 - 6*s*(p0 - p1) + (2*(sin(theta*(s - 1))*(center - p1) - sin(s*theta)*(center - p3)))/sin(theta) + (s^2*(theta^2*sin(s*theta)*(center - p3) - theta^2*sin(theta*(s - 1))*(center - p1)))/sin(theta) - (4*s*(theta*cos(s*theta)*(center - p3) - theta*cos(theta*(s - 1))*(center - p1)))/sin(theta)
%% deduce p1 == p0
clear
syms p0 p1 p3 center theta

p2 = p1
s = 0

b0 = p0*(1-s) + p1*s
b1 = p1*(1-s) + p2*s
b2 = center + (sin((1-s)*theta).*(p2 - center) + sin(s*theta).*(p3 - center)) / sin(theta)

c0 = b0*(1-s) + b1*s
c1 = b1*(1-s) + b2*s

d0 = collect(c0*(1-s) + c1*s,s)

dd0 = diff(d0,s)

%% 
c = [1,0,0]';
radius = 0.8;
p3 = c + radius*[1,0,0]';
p2 = c + radius*[cosd(30),sind(30),0]';

r1 = p3 - c;
r2 = p2 - c;

theta = acos(r1'*r2/radius/radius);

s = 0:0.001:1;

circle = c + (sin((1-s)*theta).*r1 + sin(s*theta).*r2) / sin(theta)

p_ = c + radius*[cos(s*pi)',sin(s*pi)',zeros(length(s),1)]';


hold on
plot3(circle(1,:),circle(2,:),circle(3,:))
plot3(p_(1,:),p_(2,:),p_(3,:),'--')
axis equal



%%

theta = 0:0.01:pi;

p = [...
2*cos(theta) - (theta - theta.*cos(theta).^2)./sin(theta) - 2;
                           2*sin(theta) + theta.*cos(theta)
]

r = [cos(theta);sin(theta)];
r(:,10)'*p(:,10)


