function [p,v,a,j] = s_s_curve_multiple(t, p0, v0, pos, vend, vel, acc, jerk, T, Ta, Tb, mode)
% 计算当前点位所需的最大最小时间
%
% pa     : init pos
% va     : init vel
% pb     : end  pos
% max_vb : max end vel
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period
%
% Tmax：开始时尽可能快的减速，若减速到0，则为inf，否则以到达pb的时间为准
% Tmin：开始时尽可能快的加速，直到速度最大，之后保持最大速度到终点

[m,n] = size(pos);
p = zeros(length(t), n);
for j=1:length(t)
    % 找到当前位于的区间
    t_left = t(j);
    for k=1:m
        if(t_left < T(k))
            i = k;
            break;
        end
        t_left = t_left - T(k);
    end
    for r=1:n
        if(i==1)
            p(j,r) = s_s_curve( ...
                t(j), ...
                p0(r), ...
                v0(r), ...
                pos(i,r), ...
                vend(i,r), ...
                vel(i,r), ...
                acc(i,r), ...
                jerk(i,r), ...
                T(i), ...
                Ta(i,r), ...
                Tb(i,r), ...
                mode(i,r));
        else
            p(j,r) = s_s_curve( ...
                t(j)-sum(T(1:i-1)), ...
                pos (i-1,r), ...
                vend(i-1,r), ...
                pos (i,r), ...
                vend(i,r), ...
                vel (i,r), ...
                acc (i,r), ...
                jerk(i,r), ...
                T   (i), ...
                Ta  (i,r), ...
                Tb  (i,r), ...
                mode(i,r));
        end
    end
end
end

