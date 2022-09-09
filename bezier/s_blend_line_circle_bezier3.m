function [p,dp,ddp] = s_blend_line_circle_bezier3(p0, p1, center, axis, theta, s)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
p0 = p0(:);
p1 = p1(:);
center = center(:);
axis = axis(:);

rx = p1 - center;
ry = cross(axis, rx);

% line part
line_part = (p1 - p0).*s.^3 + (3*p0 - 3*p1).*s.^2 + (3*p1 - 3*p0).*s + p0;

% circle part
circle_part = zeros(3,length(s));
for i=1:length(s)
    k = s(i);
    circle_part(:,i) = center + sin(k*k*k*theta).*ry + cos(k*k*k*theta).*rx;
end

% total sum
p  = line_part + circle_part - p1;
dp = (3*p1 - 3*p0).*s.^2 + (6*p0 - 6*p1).*s - 3*p0 + 3*p1 + ...        % line part
     3*ry.*s.^2.*theta.*cos(s.^3.*theta) - 3*rx.*s.^2.*theta.*sin(s.^3.*theta); % circle part

s3t = s.^3.*theta;
co = cos(s3t);
si = sin(s3t);
k1 = 6*s.*theta;
k2 = 9*s.^4.*theta.^2;

ddp= (6*p1 - 6*p0).*s + 6.*p0 - 6.*p1 +...
     ry.*k1.*co - ry.*k2.*si - rx.*k2.*co - rx.*k1.*si;

end