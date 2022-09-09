function [p,dp,ddp] = s_blend_circle_circle_bezier3(p1, c1, ax1, theta1, c2, ax2, theta2, s)
%UNTITLED 此处提供此函数的摘要
%   p1   : 两个圆的交点
%   c1   : circle1 的圆心坐标
%   c2   : circle2 的圆心坐标
% theta1 : 第一个圆弧的转动角度
% theta2 : 第二个圆弧的转动角度

p1 = p1(:);
c1 = c1(:);
ax1 = ax1(:);
c2 = c2(:);
ax2 = ax2(:);

% 计算p0
rx1 = p1 - c1;
ry1 = cross(ax1, rx1);

rx2 = p1 - c2;
ry2 = cross(ax2, rx2);

circle_1 = zeros(3,length(s));
circle_2 = zeros(3,length(s));
for i=1:length(s)
    k = s(i);
    circle_1(:,i) = c1 + sin((k-1)*(k-1)*(k-1)*theta1).*ry1 + cos((k-1)*(k-1)*(k-1)*theta1).*rx1;
    circle_2(:,i) = c2 + sin(k*k*k*theta2).*ry2 + cos(k*k*k*theta2).*rx2;
end

% 计算下文所需的一些临时变量
s3t1 = (s-1).^3.*theta1;
co1 = cos(s3t1);
si1 = sin(s3t1);
j1 = 6*(s-1).*theta1;
k1 = 9*(s-1).^4.*theta1.^2;

s3t2 = s.^3.*theta2;
co2 = cos(s3t2);
si2 = sin(s3t2);
j2 = 6*s.*theta2;
k2 = 9*s.^4.*theta2.^2;

% total sum
p  = circle_1 + circle_2 - p1;
dp = 3*ry1.*(s - 1).^2.*theta1.*co1 - 3.*rx1.*(s - 1).^2.*theta1.*si1 + ...
     3*ry2.*      s.^2.*theta2.*co2 - 3.*rx2.*      s.^2.*theta2.*si2;

ddp = ry1.*j1.*co1 - ry1.*k1.*si1 - rx1.*k1.*co1 - rx1.*j1.*si1 +...
      ry2.*j2.*co2 - ry2.*k2.*si2 - rx2.*k2.*co2 - rx2.*j2.*si2;
end