function [vb, v, Ta, Tb] = s_make_s_curve(pa, va, pb, max_v, a, j, T)
% 在给定过程中最大速度v，加速度a，跃度j，时间长度T的情况下
% 自动计算在最小的末端速度vb的情况下，整个过程的s曲线
%
% t      : current time
% pa     : init pos
% va     : init vel
% pb     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period

cons = eps * 10;

pt = pb - pa;

% 计算vb_max & vb_min
vb_max = min(max_v, s_acc_vend(va, a, j, T));
vb_min = max(0    , s_acc_vend(va,-a,-j, T)); % vb 一定小于max_v

% 即使结束速度和中间速度为最大，也无法在指定时间内走完
if(s_curve_length(va, vb_max, max_v, a, j, T) < pt - 1e-10 && s_curve_length(va, vb_max, max_v, a, j, T) >= 0)
    s_curve_length(va, vb_max, max_v, a, j, T)
    pt
    T
    va
    vb_max
    max_v
    a
    j
%     va

    
    error('failed to proceed 1');
end

% 即使结束速度和中间速度为最小，也会超出范围
if(s_curve_length(va, vb_min, 0, a, j, T) > pt + 1e-10 || s_curve_length(va, vb_min, 0, a, j, T) < 0)
    T
    s_curve_length(va, vb_min, 0, a, j, T)
    pt
    
    error('failed to proceed 2');
end


Ta_ = s_acc_time(va,vb_max,a,j);
vb_ = s_acc_vend(vb_max,-a,-j,T-Ta_);

% Ta : 第一段加减速时间
% Tb : 第二段加减速时间
% Tc : 匀速段时间
%
% T1 : a&b段时间：正好没有匀加速段的时间（可以发生在a、b段）为：2*a/j
% T2 : a段时间：此时正好从va加速到v2，v2为b段正好没有匀加速段降到0时的速度，
% T3 : a段过程，a段正好加速到 max_v 所用时间
% T4 : a段时间：正好将 va    减速到 0
% T5 : a段时间：正好将 va    减速到 vb
% T6 : b段时间：此时正好将 max_v 减速到 0
% T7 : b段时间：此时正好将 max_v 减速到 vb
% T8 : b段时间：此时正好将 v8 减速到0，v8为a段正好没有匀加速段时的终止速度
% 

% T1 = s_acc_time(va,     0, a, j)
% T2 = 2*a/j
% T3 = s_acc_time(va, max_v, a, j)
% T4 = s_acc_time(va,     0, a, j)
% T5 = s_acc_time(va,    vb, a, j)
% T6 = s_acc_time(max_v,  0, a, j)
% T7 = s_acc_time(max_v, vb, a, j)
% T8 = s_acc_time(v8   ,  0, a, j)

% 控制距离：
% ------------------ l1 --------------------- 
% 全力减速
Ta = min(T,s_acc_time(va,     0, a, j));
Tb = 0;
Tc = T-Ta;
vb = s_acc_vend(va,a,j,Ta);
v  = vb;
l1 = Ta * (va + vb) / 2.0;
% ------------------ l2 --------------------- 
% 先加速到最大加速度处，再全力减速
Ta = min([T,2*a/j,s_acc_time(va, max_v, a, j)]);
v  = s_acc_vend(va,         a, j,Ta);
Tb = min(T-Ta, s_acc_time(v, 0, a, j));
Tc = 0;
vb = s_acc_vend(va + a^2/j,-a,-j,Tb);
l2 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0;
% ------------------ l3 --------------------- 
% 先加速到最大速度处，再全力减速
% Ta = min(T,s_acc_time(va, max_v, a, j))
% v  = s_acc_vend(va, a, j, Ta)
% Tb = min(T-Ta, s_acc_time(v, 0, a, j))
% Tc = 0
% vb = s_acc_vend(v ,-a,-j,Tb)
% l3 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0
% ------------------ l4 --------------------- 
% 先加速到最大速度处，在末端减速到0的情况下，尽可能长的匀速
% Ta = s_acc_time(va, max_v, a, j)
% Tb = min(T-Ta, s_acc_time(max_v, 0, a, j))
% Tc = T-Ta-Tb
% v  = max_v
% vb = 0
% l4 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0 + Tc * max_v
% ------------------ l5 --------------------- 
% 先加速到最大速度处，之后一直匀速
% Ta = s_acc_time(va, max_v, a, j)
% Tb = 0
% Tc = T-Ta-Tb
% v  = max_v
% vb = max_v
% l5 = Ta * (va + v) / 2.0 + Tc * max_v

% CASE 1: 0  <  pt < p1 
%   无法行进
% CASE 2: p1 <= pt < p2 
%   p2为潜在可能的最长行进距离


% CASE 1: 0  <  pt < p1 
%   无法行进
% CASE 2: p1 <= pt < p2
%   v  = va + j*Ta*Ta/4
%   la = j/8*Ta^3 + va*Ta
%   【SUBCASE 1】: 可减速到0, 且无匀减速段
%     Tb = sqrt(4*va/j + Ta*Ta)
%     lb = Tb*v/2
%     lc = T - Ta - Tb
%     【condition】: 
%     Ta + Tb <= T
%     Tb <= 2*a/j
%     =>:
%     Ta <= (T^2-4*va/j)/(2*T)
%     Ta <= (2*sqrt(a^2 - j*va))/j
%     =>:
%     3*Ta^2 + 2*T*Ta + (4*va)/j - T^2 <= 0
%     Ta^2                             <= a^2/j^2-va/j
%     =>:
%     0 <= Ta <= r
%     其中:
%     r  = min(r1,r2)
%     r1 = (T^2-4*va/j)/(2*T);
%     r2 = (2*sqrt(a^2 - j*va))/j;
%     
%     
%
%     此时l上限为：
%     l11 = j/8*r^3 + va*r + sqrt(4*va/j + r*r)*(va + j*r*r/4)/2
%     
%     【result】:
%
%     带入方程(pt - la)^2 - lb^2 = 0
%     可得：
%     k6*Ta^6 +k4*Ta^4 +k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
%     其中：
%     k6 = 3*j*j/64;
%     k4 = 5*j*va/16;
%     k3 = j*pt/4;
%     k2 = va*va/2;
%     k1 = 2*pt*va;
%     k0 = va^3/j - pt^2;
%     
%     基于newton法求解上述方程
%   【SUBCASE 2】: 可减速到0, 且有匀减速段
%     Tb = v / a + a/j;
%     lb = Tb*v/2
%
%     【condition】: 
%     Ta + Tb <= T
%     Tb > 2*a/j
%     =>
%     j/4/a*Ta^2 + Ta - T + a/j + va/a <= 0
%     j/4/a*Ta^2 > a/j - va/a
%
%     【result】:
%     带入方程la + lb = pt
%     可得：
%     k4*Ta^4 +k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
%     其中：
%     k4 = j^2/(32*a)
%     k3 = j/8
%     k2 = ((j*(a/j + va/a))/8 + (j*va)/(8*a))
%     k1 = va
%     k0 = (va*(a/j + va/a))/2 - pt
%
%   【SUBCASE 3】: 不可减速到0, 且无匀减速段
%     Tb = T-Ta
%     lb = Tb*v/2
%
%     带入方程la + lb = pt
%     可得：
%     k2*Ta^2 +k1*Ta + k0
%     其中：
%     k2 = T*j/8
%     k1 = va/2
%     k0 = T*va/2 - pt
%
%   【SUBCASE 4】: 不可减速到0, 且有匀减速段
%     Tb = T-Ta
%     lb = Tb*v/2
%
%     带入方程la + lb = pt
%     可得：
%     k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
%     其中：
%     k3 = - j/8
%     k2 = T*j/4 - 5/2
%     k1 = 5*T - 25/(2*j)
%     k0 = T*(2*va - 5*T + 25/j)/2 - pt

Ta_max = s_acc_time(va, max_v, a, j);
Tb_max = s_acc_time(0, max_v, a, j);
pacc_max = Ta_max*(va+max_v)/2;
pdec_max = Tb_max*(0+max_v)/2;

if(a^2 - j*va >= 0)
    r1  = (T^2-4*va/j)/(2*T);
    r2  = (2*sqrt(a^2 - j*va))/j;
    r   = min(r1,r2);
    l11 = j/8*r^3 + va*r + sqrt(4*va/j + r*r)*(va + j*r*r/4)/2;
else
    l11 = -1;
end

% 二分法找vb，尽可能让结束时的速度最小
if(pt < l1)
elseif(pt < l2 && pt < l11)
    k4 = pt^2 - va^3/j;
    k3 = -2*pt*va;
    k2 = va^2/4;
    k1 = -(j*pt)/4;
    k0 = (j*va)/16;

    Ta  = r/2;
    dV  = 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta + k1;
    dTa = (k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k1*Ta + k0)/dV;
    Ta_next  = Ta - dTa;
    while(abs(Ta_next - Ta) > 1e-14)
        Ta  = Ta_next;
        dV  = 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta + k1;
        dTa = (k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k1*Ta + k0)/dV;
        Ta_next  = Ta - dTa;
    end
    Ta = Ta_next;
    vb = 0;
% elseif(pt < l2 && pt < l11)
%     k6 = 3*j*j/64;
%     k4 = 5*j*va/16;
%     k3 = j*pt/4;
%     k2 = va*va/2;
%     k1 = 2*pt*va;
%     k0 = va^3/j - pt^2;
% 
%     Ta  = r/2;
%     dV  = 6*k6*Ta^5 + 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta + k1;
%     dTa = (k6*Ta^6 + k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k1*Ta + k0)/dV;
%     Ta_next  = Ta - dTa;
%     while(abs(Ta_next - Ta) > 1e-14)
%         Ta  = Ta_next;
%         dV  = 6*k6*Ta^5 + 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta + k1;
%         dTa = (k6*Ta^6 + k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k1*Ta + k0)/dV;
%         Ta_next  = Ta - dTa;
%     end
%     Ta = Ta_next;
%     vb = 0;
elseif(s_curve_length(va, vb_min, max_v, a, j, T) > pt)
    % 近似 CASE 1
    vb = vb_min;
% elseif(Ta_ * (va + vb_max) + vb_max * (T - Ta_) < pt)
%     vb = vb_max;
elseif(Ta_max + Tb_max < T ...
    && s_curve_length(va, 0, max_v, a, j, T) < pt)
    
    if(pacc_max - a^3/j^2 + max_v*(T - Ta_max) < pt)
        Tb = ((pacc_max + max_v*(T - Ta_max) - pt)*8/j)^(1/3);
        vb = max_v - j*Tb*Tb/4;
    else
        B = -a/j;
        C = -(2*pacc_max - 2*pt + 2*max_v*(T - Ta_max))/a;
        Tb = (-B+sqrt(B^2-4*C))/2;
        vb = max_v - Tb*a + a^2/j;
    end
else
    vb_failed  = vb_min;
    vb_success = vb_max;
    while(abs(vb_success - vb_failed)>cons)
        v1_next = (vb_success + vb_failed) / 2;
        if(s_curve_length(va, v1_next, max_v, a, j, T) > pt)
            vb_success = v1_next;
        else
            vb_failed  = v1_next;
        end
    end
    vb = vb_success;
end




% 让整个过程的速度v尽可能小
% 二分法找v
v_upper = max_v;
v_below = 0;
while(abs(v_upper - v_below)>cons)
    v_next = (v_upper + v_below) / 2;
    
    [preal, vreal] = s_curve_length(va, vb, v_next, a, j, T);

    if(preal > pt)
        v_upper = v_next;
    else
        v_below = v_next;
    end
end
v=vreal;

% finally gives result
if(v >= max(va,vb) || v <= min(va,vb))
    Ta = s_acc_time(va,v,a,j);
    Tb = s_acc_time(vb,v,a,j);
else
    total_even_time = T - s_acc_time(va,vb,a,j);
    Ta = (1-abs(va-v)/abs(va-vb)) * total_even_time;
    Tb = (1-abs(vb-v)/abs(va-vb)) * total_even_time;
end

end





