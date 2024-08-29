function [T, Ta, Tb, vb, vc, mode] = s_make_s_curve_multiple(p0, v0, pos, vb_max, vc_max, acc, jerk)
% 计算当前点位所需的最大最小时间
%
% pa        : init pos
% va        : init vel
% pb        : end  pos
% max_vend  : max end vel
% max_vel   : max vel  during period
% acc       : max acc  during period
% jerk      : max jerk during period
%
% T         : time
% Ta        : time acc
% Tb        : time dec
% real_vend : real velocity at end of each path 
% real_vel  : real max velocity of each path

[m,n] = size(pos);

Tmax      = zeros(m,n);
Tmin      = zeros(m,n);
T         = zeros(m,1);
vb = zeros(m,n);
vc  = zeros(m,n);
Ta        = zeros(m,n);
Tb        = zeros(m,n);
mode      = zeros(m,n);

for i = 1:m
%     i
    for j = 1:n
%         j
        if(i == 1)
            pa = p0(1,j);
            va = v0(1,j);
        else
            pa = pos(i-1,j);
            va = v0(j);
        end
        [Tmax(i,j), Tmin(i,j)] = s_s_curve_Tmax_Tmin( ...
            pa, ...
            va, ...
            pos(i,j), ...
            vb_max(i,j), ...
            vc_max(i,j), ...
            acc(i,j), ...
            jerk(i,j));
    end
%     pa
%     va
%     pos
%     max_vend
%     max_vel
%     acc
%     jerk
%     Tmax
%     Tmin
    
    if(i==7)
%         i
    end
    
    Tmax_all = min(Tmax(i,:));
    Tmin_all = max(Tmin(i,:));
    
    % 二分法 search max T
    if(Tmax_all == inf)
        T_upper = 1000;
    else
        T_upper = Tmax_all;
    end
    T_below = Tmin_all;
    
    diff      = abs(T_upper - T_below);
    diff_last = diff * 2;
    while(diff < diff_last)
        diff_last = diff;

        T_next = (T_upper + T_below)/2;
        
        if(s_test_curve_slow(i, p0, v0, pos, vb_max, vc_max, acc, jerk, T_next))
            T_upper = T_next;
        else
            T_below = T_next;
        end

        diff = abs(T_upper - T_below);
    end
    T(i) = T_upper;
    
    if(i==19)
        i;
    end
    for j=1:n
        [vb(i,j),vc(i,j),Ta(i,j),Tb(i,j),mode(i,j)]=s_make_s_curve( ...
            p0(j), ...
            v0(j), ...
            pos(i,j), ...
            vc_max(i,j), ...
            acc(i,j), ...
            jerk(i,j), ...
            T(i));
%         p0
%         v0
%         pos
%         max_vel
%         acc
%         jerk
%         T

    end
    p0 = pos(i,:);
    v0 = vb(i,:);
end

end

