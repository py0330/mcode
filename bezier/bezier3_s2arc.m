function [arc, darc] = bezier3_s2arc(p0,p1,p2,s)
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


% x2 = sin(theta) + 1;
% y2 = cos(theta);

% p0 = [0,0,0]';
% p1 = [1,0,0]';
% p2 = [x2,y2,0]';


a = 12-6*cos(theta/2);
b = -a;
c = 3;

d = (a-6)/pi;
e = (a-6)/pi*0.5;
f = -(a-6)/pi*0.5;

darc = dis*(a*s.*s + b*s + c + d*sin(s*pi) + e*cos(s*2*pi) + f);
arc  = dis*(1/3*a*s.*s.*s + 1/2*b*s.*s - d/pi*cos(s*pi) + e/2/pi*sin(s*2*pi) + (f+c)*s + d/pi);
end








