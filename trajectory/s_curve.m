function [p, v1] = s_curve(t, p0, v0, p1, v, a, j, T)
% 在给定过程中最大速度v，加速度a，跃度j，时间长度T的情况下
% 自动计算在最小的末端速度v1的情况下，整个过程的s曲线
%
% t      : current time
% p0     : init pos
% v0     : init vel
% p1     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period

cons = eps * 10;

% si 如下：si = sign(p1 - T*v0)     |   若（p1 - T*v0 > 0）则 v1 > v0
pt = p1 - p0;

% 计算v1_max & v1_min
v1_max = min(v, s_acc_vend(v0,a,j,T));
v1_min = max(0, s_acc_vend(v0,-a,-j,T));

% 即使结束速度和中间速度为最大，也无法在指定时间内走完
if(s_curve_length(v0, v1_max, v, a, j, T) < pt)
    error('failed to proceed 1');
end

% 即使结束速度和中间速度为最小，也会超出范围
if(s_curve_length(v0, v1_min, 0, a, j, T) > pt)
    error('failed to proceed 2');
end

% 二分法找v1，尽可能让结束时的速度最小
if(s_curve_length(v0, v1_min, v, a, j, T) > pt)
    v1 = v1_min;
else
    v1_failed  = v1_min;
    v1_success = v1_max;
    while(abs(v1_success - v1_failed)>cons)
        v1_next = (v1_success + v1_failed) / 2;
        if(s_curve_length(v0, v1_next, v, a, j, T) > pt)
            v1_success = v1_next;
        else
            v1_failed  = v1_next;
        end
    end
    v1 = v1_success;
end

% 二分法找v，让整个过程的速度尽可能小
v_upper = v;
v_bottom = 0;
while(abs(v_upper - v_bottom)>cons)
    v_next = (v_upper + v_bottom) / 2;
    if(s_curve_length(v0, v1, v_next, a, j, T) > pt)
        v_upper = v_next;
    else
        v_bottom  = v_next;
    end
end
v=(v_upper+v_bottom)/2;

% finally calculate p
Ta = s_acc_time(v0,v,a,j);
Tb = s_acc_time(v1,v,a,j);

if(t < Ta)
    si = sign(v - v0);

    if(Ta >= 2*a/j)
        if(t < a/j)
            p = p0 + v0 * t + si/6*j*t^3;
        elseif(t < (Ta - a/j))
            p = p0 + v0 * a / j + si/6 * a^3 / j^2 + (v0 + si/2*a^2/j) * (t - a/j)+ si*a/2*(t-a/j)^2; 
        else
            p = p0 + (v0+v)/2*Ta-(v*(Ta-t) - si/6*j*(Ta-t)^3);
        end
    else
        if(t < Ta / 2)
            p = p0 + v0 * t + si/6*j*t^3;
        else
            p = p0 + (v0+v)/2*Ta - (v * (Ta-t) - si/6*j*(Ta-t)^3);
        end
    end
elseif(t < T-Tb)
    p = p0 + (v0+v)/2*Ta + v * (t - Ta);
else
    si = sign(v1 - v);
    if(Tb >= 2*a/j)
        if(T-t < a/j)
            p = p1 - v1 * (T-t) + si/6*j*(T-t)^3;
        elseif(T-t < (Tb - a/j))
            p = p1 - v1 * a / j + si/6 * a^3 / j^2 - (v1 - si/2*a^2/j) * (T - t - a/j) + si*a/2*(T - t - a/j)^2;
        else
            p = p1 - (v+v1)/2*Tb + (v*(Tb-T+t) + si/6*j*(Tb-T+t)^3);
        end
    else
        if(T-t < Tb / 2)
            p = p1 - v1 * (T-t) + si/6*j*(T-t)^3;
        else
            p = p1 - (v+v1)/2*Tb + (v*(Tb-T+t) + si/6*j*(Tb-T+t)^3);
        end
    end
end



% p = (v0 + v)/2 * T0 + (v1 + v)/2 * T1 + v * (T - T0 - T1);




end

function v = s_acc_vend(v0,a,j,T)
    if(a/j > T/2)
        v = v0 + j*T*T/4;
    else
        v = v0 + (T-a/j*2)*a + a^2/j;
    end
end

function t = s_acc_time(v0,v1,a,j)
    % CASE 1 |v - v0| >  a^2 / j, 此时加速段可以达到最大加速度：
    %        T0 = |v - v0|/a + a/j            |    |v - v0| >  a^2 / j
    % CASE 2 |v - v0| <= a^2 / j, 此时加速段（或减速段）无法达到最大加速度：
    %        T1 = 2 * sqrt( |v - v0| / j )    |    |v - v0| <= a^2 / j

    if(abs(v1-v0)>a^2 / j)
        t = abs(v1-v0) / a + a/j;
    else
        t = 2 * sqrt( abs(v1-v0) / j );
    end
end

function p = s_curve_length(v0, v1, v, a, j, T)
    cons = eps * 10;

    % 确定根据起始终止速度要求，可以满足时间T
    if(s_acc_time(v0,v1,a,j) > T)
        error('failed in curve length');
    end
    
    % 确定是否可以达到指定的v
    if(s_acc_time(v0,v,a,j) + s_acc_time(v1,v,a,j) > T)
        if(v > (v0+v1)/2)
            v_success = max(v0,v1);
            v_failed  = v;
        else
            v_success = min(v0,v1);
            v_failed  = v;
        end
        
        while(abs(v_success - v_failed)>cons)
            v_next = (v_success + v_failed)/2;
            if(s_acc_time(v0,v_next,a,j) + s_acc_time(v1,v_next,a,j) > T)
                v_failed  = v_next;
            else
                v_success = v_next;
            end
        end
        v = v_success;
    end

    T0 = s_acc_time(v0,v,a,j);
    T1 = s_acc_time(v1,v,a,j);
    
    p = (v0 + v)/2 * T0 + (v1 + v)/2 * T1 + v * (T - T0 - T1);
end
