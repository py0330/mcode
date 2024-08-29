function [dis, r] = s_is_in_vavg_boundage(v0,v1,a,dt,vavg)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

v0 = v0(:);
v1 = v1(:);
vavg = vavg(:);

t1 = norm(v1 - v0) / a;
r = (dt - t1)/2*a;

if(abs(dt - t1) < 1e-10)
    dis = norm(vavg - (v0 + v1)/2);
    return;
end

vavg_left = (vavg*dt - (v0 + v1)/2*t1)/(dt-t1);
% vavg=vavg_left;

a = vavg_left - v0;
b = vavg_left - v1;
c = v1 - v0;

if(norm(c) < 1e-10 || a'*c < 0)
    dis = norm(a);
elseif(-b'*c <0)
    dis = norm(b);
else
    % 点到直线的距离计算公式
    dis = norm(cross(vavg_left - v0, v1-v0))/norm(v1-v0);
end




end