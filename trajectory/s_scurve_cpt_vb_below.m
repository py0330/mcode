function [vb_below] = s_scurve_cpt_vb_below(pa, va, pb, vc_max, vb_max, a, j, T)
%
% 1. 确定必然可以达到的 vb_max vb_min vc_max vc_min
%
% 2. 确定 vb_min 对应的最大可能的vc：vc_l（vb_max 对应 vc_s）
% 
% 3. 确定 vc_max 对应的最小可能vb：vb_s（vc_min 对应 vb_l）
%    
%    vc        vb           Ta         Tc       Tb         requires
% 1 vc_max   vb_max       T_va2vcmax T-Ta-Tb T_vcmax2vbmax none
% 2 vc_max   vc_max-a*a/j T_va2vcmax T-Ta-Tb 2*a/j         T, vb
% 3 vc_max   vb_s         T_va2vcmax T-Ta-Tb T_vcmax2vbs   none    
% 4 va/-Ta   vc-a*a/j     T-Tb       0       2*a/j         T > 4*a/j
% 5 va+a*a/j vc\-Tb       2*a/j      0       T-Ta          T > 4*a/j
% 6 vc_l     vb_min        Ta        Tc         Tb         none
%
% l1 -> l3 -> l6 依次减小
%
% 上式先计算 l1 l3 l6
% l2 若 vc_max-a*a/j > vb_max, 则 l2 = l1
%    若 vc_max-a*a/j < vb_min, 则 l2 = l3
%    否则 按照公式计算
%
% l4 按照公式计算
%
% l5 按照公式计算
% 
% l1 <= pt      时：vb_upper = vb_max
% l2 <= pt < l1 时：有匀速段，Tb段无匀加速，计算Tb
% l3 <= pt < l2 时：有匀速段，Tb段有匀加速，计算Tb
%   T > 4*a/j 时：l4 > l5
%     l4 <= pt < l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l5 <= pt < l4 时：无匀速段，Ta段有匀加速，Tb段有匀加速
%     l6 <= pt < l5 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 4*a/j && T > 2*a/j 时：l5 > l4
%     l5 <= pt < l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l4 <= pt < l5 时：无匀速段，Ta段无匀加速，Tb段无匀加速
%     l6 <= pt < l4 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 2*a/j 时：
%     l6 <= pt < l3 时：无匀速段，Ta段无匀加速，Tb段无匀加速
% pt < l6 时：vb_upper = vb_min
cons = 1e-10;
Z1 = a^2/j;
Z2 = T^2*j;
pt = pb - pa;

vb_min = 0;

% STEP 1: 确定必然可以达到的 vb_max vb_min vc_max vc_min
vb_max = min(vb_max, s_acc_vend(va,a,j,T));
vb_min = max(vb_min, s_acc_vend(va,-a,-j,T));
vc_max = min(vc_max, s_cpt_vc_upper_by_va_vb_T(va,vb_max,T,a,j));

T_va_to_vcmax = s_acc_time(va,vc_max,a,j);
T_vcmax_to_vbmax = s_acc_time(vc_max,vb_max,a,j);

% STEP 2: 确定 vb_min 对应的最大可能的vc：vc_l
vc_l = min(vc_max, s_cpt_vc_upper_by_va_vb_T(va,vb_min,T,a,j));

% STEP 3: 确定 vc_max 对应的最小可能vb：vb_s
vb_s = max(vb_min, s_acc_vend(vc_max,-a,-j,T-T_va_to_vcmax));

% CASE 1 ------------------ pt > l1
vc = vc_max;
vb = vb_max;
Ta = T_va_to_vcmax;
Tb = T_vcmax_to_vbmax;
Tc = T - Ta - Tb;
l1 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
if(pt >= l1 + 1e-10)
    vb_below = inf;
    return;
end
if(pt >= l1)
    vb_below = vb_max;
    return;
end

% CASE 2 ------------------ l3 < pt < l1
vc = vc_max;
vb = vb_s;
Ta = T_va_to_vcmax;
Tb = s_acc_time(vc,vb,a,j);
Tc = T - Ta - Tb;
l3 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
if(pt >= l3)
    if(vc_max-a*a/j > vb_max)
        l2 = l1;
    elseif(vc_max-a*a/j < vb_s)
        l2 = l3;
    else
        vc = vc_max;
        vb = vc_max-a*a/j;
        Ta = T_va_to_vcmax;
        Tb = s_acc_time(vc, vb, a,j);
        Tc = T - Ta - Tb;
        l2 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
    end
    if(pt >= l2)
        % clear
        % syms la vc_max pt Ta Tb T
        % vb = vc_max - j*Tb*Tb/4;
        % l = la + (T - Ta - Tb)*vc_max + (vc_max + vb)/2*Tb
        % solve(l==pt, Tb)
        Ta = T_va_to_vcmax;
        la = Ta*(va+vc_max)/2;
        Tb = ((la + vc_max*(T - Ta) - pt)*8/j)^(1/3);
        vb = vc_max - j*Tb*Tb/4;
    
        vb_below = vb;
        return;
    else
        Ta = T_va_to_vcmax;
        la = Ta*(va+vc_max)/2;
        B = -a/j;
        C = -(2*la - 2*pt + 2*vc_max*(T - Ta))/a;
    
        Tb = (-B+sqrt(max(B^2-4*C,0)))/2;
        vb = vc_max - Tb*a + Z1;
    
        vb_below = vb;
        return;
    end
end

% CASE 3 ------------------ l6 < pt < l3
vc = vc_l;
vb = vb_min;
Ta = s_acc_time(va, vc, a,j);
Tb = s_acc_time(vc, vb, a,j);
Tc = T - Ta - Tb;
l6 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
if(pt >= l6)
    Tb = 2*a/j;
    Ta = T-Tb;
    vc = s_acc_vend(va,a,j,Ta);
    vb = s_acc_vend(vc,-a,-j,Tb);
    l4 = (va+vc)*Ta/2 + (vc+vb)*Tb/2;

    Ta = 2*a/j;
    Tb = T-Ta;
    vc = s_acc_vend(va,a,j,Ta);
    vb = s_acc_vend(vc,-a,-j,Tb);
    l5 = (va+vc)*Ta/2 + (vc+vb)*Tb/2;

    % Ta无匀加速, Tb无匀加速
    if((T <= 2*a/j) || (2*a/j < T && T <= 4*a/j && l4 <= pt && pt < l5))
        % clear
        % syms va j Ta a T pt
        % vc = va + j*Ta*Ta/4
        % la = Ta*(va+vc)/2
        % Tb = T-Ta
        % vb = vc - j*Tb*Tb/4
        % lb = Tb*(vc+vb)/2
        % collect(la+lb-pt,Ta)
        % 【result】:
        % 带入方程la + lb = pt
        % 可得：
        % k2*Ta^2 +k1*Ta + k0
        % 其中：
        % k2 = T*j/8
        % k1 = -3*T^2*j/8
        % k0 = pt - (T*(2*va - (T^2*j)/4))/2
        %
        % Ta = (-k1+sqrt(k1*k1-4*k0*k2))/2/k2
        k2 = T*j/8;
        k1 = -3*T^2*j/8;
        k0 = pt - (T*(2*va - T^2*j/4))/2;
    
        % 选根
        % 其极值为 (k1)/(2*k2) = (3*T)/2
        % 因此需选其较小的根
        Ta = (-k1-sqrt(k1*k1-4*k0*k2))/2/k2;
    
        vc = va + j*Ta*Ta/4;
        Tb = T-Ta;
        vb_below = vc - j*Tb*Tb/4;
        vb = vb_below;
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < -cons || vc > vc_max + cons || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.1')
        end
        return;
    end
    
    % Ta有匀加速, Tb无匀加速
    if((4*a/j < T && l4 <= pt && pt < l3) ||...
            (2*a/j < T && T <= 4*a/j && l5 <= pt && pt < l3))
        % clear
        % syms va T a j pt Ta
        % vc = va + Ta*a - a^2/j
        % la = Ta*(va + vc)/2
        % Tb = T-Ta
        % vb = vc - j*Tb*Tb/4;
        % lb = Tb*(vc + vb)/2
        % collect(la+lb-pt,Ta)
        %
        %   【result】:
        %   带入方程la + lb = pt
        %   可得：
        %   k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
        %   其中：
        %   k3 = j/8
        %   k2 = (- a/2 - (3*T*j)/8)
        %   k1 = ((T^2*j)/8 + a^2/(2*j) + (T*(2*a + (T*j)/2))/2)
        %   k0 = - pt - (T*((T^2*j)/4 - 2*va + (2*a^2)/j))/2
    
        %
        %   【condition】:
        %   Tb > 2*a/j
        %   => Ta < T - 2*a/j
        %   于是：
        %   0 <= Ta <= min(T, 2*a/j, T - 2*a/j)
        %
        %   计算Ta = min(T, 2*a/j, T - 2*a/j)
        %   l3 = la + lb
        k3 = j/8;
        k2 = (- a/2 - (3*T*j)/8);
        k1 = ((T^2*j)/8 + a^2/(2*j) + (T*(2*a + (T*j)/2))/2);
        k0 = - pt - (T*((T^2*j)/4 - 2*va + (2*a^2)/j))/2;
        
        % 选根
        % syms f(Ta) g(Ta)
        % f(Ta) = k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
        % g(Ta) = diff(f,Ta)
        % solve(g,Ta)
        %
        % 可得其极值：
        % r1 = 4/3*T - 2/3*a/j
        % r2 = -(2*a)/j
        %
        % 由于b段可达最大加速度，因此必有 T >= 2a/j
        % 于是
        % r1 >= 4/3*T - 1/3*T = T
        % 因此其上下界为 [0,T]
        Ta = newton_raphson_binary_search(@(x)(k3*x*x*x+k2*x*x+k1*x+k0),0,T,10*eps);
    
        vc  = va + Ta*a - Z1;
        Tb = T - Ta;
        vb = s_acc_vend(vc,-a,-j,Tb);
        vb_below = vb;
    
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < -cons || vc > vc_max + cons || abs(l-pt) > max(pt,1) * cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.3')
        end
    
        return;
    
    end
    
    % Ta无匀加速, Tb有匀加速
    if((4*a/j < T && l6 <= pt && pt < l5) ||...
            (2*a/j < T && T <= 4*a/j && l6 <= pt && pt < l4))
        %   clear
        %   syms va j Ta a T
        %   vc = va + j*Ta*Ta/4
        %   la = Ta*(vc+va)/2
        %   Tb = T-Ta
        %   vb = vc - Tb*a + a^2/j;
        %   lb = Tb*(vc+vb)/2
        %   l  = la + lb
        %
        %   【result】:
        %   带入方程la + lb = pt
        %   可得：
        %   k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
        %   其中：
        %   k3 = - j/8
        %   k2 = (T*j)/4 - a/2
        %   k1 = - a^2/(2*j) + T*a
        %   k0 = (T*(a^2/j - T*a + 2*va))/2 - pt
        %
        %   【condition】:
        %   Tb > 2*a/j
        %   => Ta < T - 2*a/j
        %   于是：
        %   0 <= Ta <= min(T, 2*a/j, T - 2*a/j)
        %
        %   计算Ta = min(T, 2*a/j, T - 2*a/j)
        %   l3 = la + lb
        k3 = - j/8;
        k2 = (T*j)/4 - a/2;
        k1 = - Z1/2 + T*a;
        k0 = (T*(Z1 - T*a + 2*va))/2 - pt;
    
        % 选根
        % syms f(Ta) g(Ta)
        % f(Ta) = k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
        % g(Ta) = diff(f,Ta)
        % solve(g,Ta)
        %
        % 可得其极值：
        % r1 = 4/3*T - 2/3*a/j
        % r2 = -(2*a)/j
        %
        % 由于b段可达最大加速度，因此必有 T >= 2a/j
        % 于是
        % r1 >= 4/3*T - 1/3*T = T
        % 因此其上下界为 [0,T]
        Ta = newton_raphson_binary_search(@(x)(k3*x*x*x+k2*x*x+k1*x+k0),0,T,10*eps);
    
        vc = va + j*Ta*Ta/4;
        Tb = T - Ta;
        vb_below = s_acc_vend(vc,-a,-j,Tb);
        vb = vb_below;
    
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < -cons || vc > vc_max + cons || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.2')
        end
    
        return;
    
    
    end
    
    % Ta有匀加速, Tb有匀加速
    if((4*a/j < T && l5 <= pt && pt < l4))
        % syms va T a j pt Ta
        % vc = va + Ta*a - a^2/j
        % la = Ta*(va + vc)/2
        % Tb = T-Ta
        % vb = vc- Tb*a + a^2/j;
        % lb = Tb*(vc+vb)/2
        %
        % 根据 la + lb = pt，有：
        % collect(la+lb-pt,Ta)
        % - a*Ta^2 + (2*T*a + a^2/j)*Ta - pt - (T*(T*a - 2*va + (3*a^2)/j))/2
        %
        % 可得方程系数
        k2 = - a;
        k1 = 2*T*a;
        k0 = - pt - (T*(T*a - 2*va + a^2/j))/2;
    
        % 选根
        % syms f(Ta) g(Ta)
        % f(Ta) = k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
        % g(Ta) = diff(f,Ta)
        % solve(g,Ta)
        %
        % 得到:
        %
        % r1 = T + (2*a)/(3*j)
        % r2 = T + (2*a)/j
        %
        % 因为必有 2*a/j <= Ta <= T, 因此其上下界为：
        % T_below = 2*a/j
        % T_upper = T
    
        Ta = (-k1+sqrt(k1*k1-4*k2*k0))/(2*k2);
    
        vc  = va + Ta*a - Z1;
        Tb = max(T - Ta, 0);
        vb = s_acc_vend(vc,-a,-j,Tb);
        vb_below = vb;
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < -cons || vc > vc_max + cons || abs(l-pt) > cons || abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.4')
        end
    
        return;
    end
end

% CASE 4 ------------------ pt < l6
if(pt < l6)
    vb_below = vb_min;
    return;
end



error('condition check failed');
end