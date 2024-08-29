function [vavg_max,vavg_min] = s_v_avg_boundage(v0, v1, a, dt)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

% 计算给定起始与截止速度 v0 v1、区间最大加速度 a、时间周期 dt，
% 计算整个周期内的最大最小平均速度

va = min(v0, v1);
vb = max(v0, v1);

t1 = (vb - va)/a;

vavg_max = (va + vb) * t1 + (dt - t1) * (a * (dt - t1) / 2 + vb + vb) / 2;
vavg_min = (va + vb) * t1 + (dt - t1) * (-a * (dt - t1) / 2 + va + va) / 2;

end