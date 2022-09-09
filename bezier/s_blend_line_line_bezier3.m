function [p,dp,ddp] = s_blend_line_line_bezier3(p0, p1, p2, s)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
p0 = p0(:);
p1 = p1(:);
p2 = p2(:);

% 使用3阶bezier曲线，控制点为：[p0,p1,p1,p2]

p   =   (p2 - p0).*s.^3 - 3*(p1 - p0).*s.^2 + 3*(p1 - p0).*s + p0;
dp  = 3*(p2 - p0).*s.^2 - 6*(p1 - p0).*s    + 3*(p1 - p0);
ddp = 6*(p2 - p0).*s    - 6*(p1 - p0);

end