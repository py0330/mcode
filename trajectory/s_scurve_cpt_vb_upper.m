function [vb_upper] = s_scurve_cpt_vb_upper(pa, va, pb, vc_max, vb_max, a, j, T)
%
% 1. 确定必然可以达到的 vb_min vb_max vc_min vc_max
%
% 2. 确定 vb_max 对应的最大可能的vc：vc_l（vb_min 对应 vc_s）
% 
% 3. 确定 vc_min 对应的最小可能vb：vb_s（vc_max 对应 vb_l）
%    
%    vc        vb           Ta         Tc       Tb         requires
% 1 vc_min   vb_min       T_va2vcmin T-Ta-Tb T_vcmin2vbmin none
% 2 vc_min   vc_min-a*a/j T_va2vcmin T-Ta-Tb 2*a/j         T, vb
% 3 vc_min   vb_s         T_va2vcmin T-Ta-Tb T_vcmin2vbs   none    
% 4 va/-Ta   vc-a*a/j     T-Tb       0       2*a/j         T > 4*a/j
% 5 va+a*a/j vc\-Tb       2*a/j      0       T-Ta          T > 4*a/j
% 6 vc_l     vb_max        Ta        Tc         Tb         none
%
% l1 -> l3 -> l6 依次增加
%
% 上式先计算 l1 l3 l6
% l2 若 vc_min+a*a/j < vb_min, 则 l2 = l1
%    若 vc_min+a*a/j > vb_max, 则 l2 = l3
%    否则 按照公式计算
%
% l4 按照公式计算
%
% l5 按照公式计算
% 
% l1 <= pt      时：vb_upper = vb_min
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
% pt < l6 时：vb_upper = vb_max
cons = 1e-10;
Z1 = a^2/j;
Z2 = T^2*j;
pt = pb - pa;

vb_min = 0;
vc_min = 0;

% STEP 1: 确定必然可以达到的 vb_min vb_max vc_min vc_max
vb_min = max(vb_min, s_acc_vend(va,-a,-j,T));
vb_max = min(vb_max, s_acc_vend(va,a,j,T));
vc_min = max(vc_min, s_cpt_vc_below_by_va_vb_T(va,vb_min,T,a,j));

T_va_to_vcmin = s_acc_time(va,vc_min,a,j);
T_vcmin_to_vbmin = s_acc_time(vc_min,vb_min,a,j);

% STEP 2: 确定 vb_max 对应的最大可能的vc：vc_s
vc_s = max(vc_min, s_cpt_vc_below_by_va_vb_T(va,vb_max,T,a,j));

% STEP 3: 确定 vc_min 对应的最达可能vb：vb_l
vb_l = min(vb_max, s_acc_vend(vc_min,a,j,T-T_va_to_vcmin));

% CASE 1 ------------------ pt < l1
vc = vc_min;
vb = vb_min;
Ta = T_va_to_vcmin;
Tb = T_vcmin_to_vbmin;
Tc = T - Ta - Tb;
l1 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;

if(pt <= l1 - cons)
    vb_upper = -1;
    return;
end
if(pt <= l1)
    vb_upper = vb_min;
    return;
end

% CASE 2 ------------------ l1 < pt < l3
vc = vc_min;
vb = vb_l;
Ta = T_va_to_vcmin;
Tb = s_acc_time(vc,vb,a,j);
Tc = T - Ta - Tb;
l3 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
if(pt <= l3)
    % 计算 l2
    if(vc_min+a*a/j < vb_min)
        l2 = l1;
    elseif(vc_min+a*a/j > vb_l)
        l2 = l3;
    else
        vc = vc_min;
        vb = vc_min-a*a/j;
        Ta = T_va_to_vcmin;
        Tb = s_acc_time(vc, vb, a,j);
        Tc = T - Ta - Tb;
        l2 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
    end

    if(pt <= l2)
        % clear
        % syms va j Ta a T pt Tb la vc vc_min
        % vb = vc_min + j*Tb*Tb/4
        %
        % l = Tb*(vb+vc_min)/2 + la + (T-Ta-Tb)*vc_min;
        % collect(l, Tb)
        % solve(l==pt,Tb)
        Ta = T_va_to_vcmin;
        la = Ta*(va+vc_min)/2;
        Tb = (-(8*(la - pt + T*vc_min - Ta*vc_min))/j)^(1/3);
        vb = vc_min + j*Tb*Tb/4;

        vb_upper = vb;
        return;
    else
        % clear
        % syms Ta la va vc_min pt a j Tb T
        % vb = vc_min + a*Tb - a*a/j
        % l = Tb*(vb+vc_min)/2  + (T-Ta-Tb)*vc_min + la;
        % collect(l, Tb)
        % solve(l==pt,Tb)
        Ta = T_va_to_vcmin;
        la = Ta*(va+vc_min)/2;
        B = -a/j;
        C = (2*la - 2*pt + 2*vc_min*(T - Ta))/a;

        Tb = (-B+sqrt(max(B^2-4*C,0)))/2;
        vb = vc_min + Tb*a - Z1;

        vb_upper = vb;
        return;
    end
end

% CASE 3 ------------------ l3 < pt < l6
vc = vc_s;
vb = vb_max;
Ta = s_acc_time(va, vc, a, j);
Tb = s_acc_time(vc, vb, a, j);
Tc = T - Ta - Tb;
l6 = (va+vc)*Ta/2 + vc*Tc + (vc+vb)*Tb/2;
if(pt <= l6)
    Tb = 2*a/j;
    Ta = T-Tb;
    vc = s_acc_vend(va,-a,-j,Ta);
    vb = s_acc_vend(vc,a,j,Tb);
    l4 = (va+vc)*Ta/2 + (vc+vb)*Tb/2;

    Ta = 2*a/j;
    Tb = T-Ta;
    vc = s_acc_vend(va,-a,-j,Ta);
    vb = s_acc_vend(vc,a,j,Tb);
    l5 = (va+vc)*Ta/2 + (vc+vb)*Tb/2;

    % Ta无匀加速, Tb无匀加速
    if((T <= 2*a/j) || (2*a/j < T && T <= 4*a/j && l4 >= pt && pt > l5))
        % syms a j T va Ta Tb pt
        % vc = va - j*Ta*Ta/4
        % la = Ta*(va+vc)/2
        % Tb = T-Ta
        % vb = vc + j*Tb*Tb/4
        % lb = Tb*(vc+vb)/2
        % collect(la + lb - pt, Ta)
        % 【result】:
        % 带入方程la + lb = pt
        % 可得：
        % k2*Ta^2 +k1*Ta + k0
        % 其中：
        % k2 = T*j/8
        % k1 = -3*T^2*j/8
        % k0 = -pt + (T*(2*va + (T^2*j)/4))/2
        %
        % Ta = (-k1+sqrt(k1*k1-4*k0*k2))/2/k2
        k2 = T*j/8;
        k1 = -3*T*T*j/8;
        k0 = -pt + (T*(2*va + (T^2*j)/4))/2;
    
        % 选根
        % 其极值为 (k1)/(2*k2) = (3*T)/2
        % 因此需选其较小的根
        Ta = (-k1-sqrt(k1*k1-4*k0*k2))/2/k2;
    
        vc = va - j*Ta*Ta/4;
        Tb = T-Ta;
        vb_upper = vc + j*Tb*Tb/4;
        vb = vb_upper;
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < vc_min - cons || vc > vc_max + cons || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.1')
        end
        return;
    end
    
    % Ta有匀加速, Tb无匀加速
    if((4*a/j < T && l4 >= pt && pt > l3) ||...
            (2*a/j < T && T <= 4*a/j && l5 >= pt && pt > l3))
        % clear
        % syms va T a j pt Ta
        % vc = va - Ta*a + a^2/j
        % la = Ta*(va + vc)/2
        % Tb = T-Ta
        % vb = vc + j*Tb*Tb/4;
        % lb = Tb*(vc + vb)/2
        % collect(la+lb-pt,Ta)
        %
        %   【result】:
        %   带入方程la + lb = pt
        %   可得：
        %   k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
        %   其中：
        %   k3 = -j/8
        %   k2 = (a/2 + (3*T*j)/8)
        %   k1 = (- (T^2*j)/8 - a^2/(2*j) - (T*(2*a + (T*j)/2))/2)
        %   k0 = - pt + (T*(2*va + (T^2*j)/4 + (2*a^2)/j))/2
    
        %
        %   【condition】:
        %   Tb > 2*a/j
        %   => Ta < T - 2*a/j
        %   于是：
        %   0 <= Ta <= min(T, 2*a/j, T - 2*a/j)
        %
        %   计算Ta = min(T, 2*a/j, T - 2*a/j)
        %   l3 = la + lb
        k3 = -j/8;
        k2 = (a/2 + (3*T*j)/8);
        k1 = (- (T^2*j)/8 - a^2/(2*j) - (T*(2*a + (T*j)/2))/2);
        k0 = - pt + (T*(2*va + (T^2*j)/4 + (2*a^2)/j))/2;
    
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
    
        vc  = va - Ta*a + a*a/j;
        Tb = T - Ta;
        vb = s_acc_vend(vc,a,j,Tb);
        vb_upper = vb;
    
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < vc_min - cons || vc > vc_max + cons || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.3')
        end
    
        return;
    
    end
    
    % Ta无匀加速, Tb有匀加速
    if((4*a/j < T && l6 >= pt && pt > l5) ||...
            (2*a/j < T && T <= 4*a/j && l6 >= pt && pt > l4))
        % clear
        % syms va j Ta a T vc pt
        % vc = va - j*Ta*Ta/4
        % la = Ta*(va+vc)/2
        % Tb = T-Ta
        % vb = vc + Tb*a - a^2/j
        % lb = Tb*(vc+vb)/2
        % l  = la + lb
        % collect(la+lb-pt,Ta)
        %
        %  【result】:
        %   带入方程la + lb = pt
        %   可得：
        %   k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
        %   其中：
        %   k3 = j/8
        %   k2 = (a/2 - (T*j)/4)
        %   k1 = (a^2/(2*j) - T*a)
        %   k0 = - pt + (T*(2*va + T*a - a^2/j))/2
     
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
        k2 = (a/2 - (T*j)/4);
        k1 = (a^2/(2*j) - T*a);
        k0 = - pt + (T*(2*va + T*a - a^2/j))/2;
    
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
    
        vc  = va - j*Ta*Ta/4;
        Tb = T - Ta;
        vb = s_acc_vend(vc,a,j,Tb);
        vb_upper = vb;
    
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < vc_min - cons || vc > vc_max + cons || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.2')
        end
    
    
        return;
    
    
    end
    
    % Ta有匀加速, Tb有匀加速
    if((4*a/j < T && l5 >= pt && pt > l4))
        % syms va j Ta a T pt
        % vc = va - Ta*a + a*a/j
        % la = Ta*(va+vc)/2
        % Tb = T-Ta
        % vb = vc + Tb*a - a*a/j
        % lb = Tb*(vc+vb)/2;
        %
        % 根据 la + lb = pt，有：
        % collect(la+lb-pt,Ta)
        %
        % 可得方程系数
        k2 = a;
        k1 = -2*T*a;
        k0 = - pt + (T*(2*va + T*a + a^2/j))/2;
    
        Ta = (-k1-sqrt(k1*k1-4*k2*k0))/(2*k2);
    
        vc  = va - Ta*a + a*a/j;
        Tb = max(T - Ta, 0);
        vb = s_acc_vend(vc,a,j,Tb);
        vb_upper = s_acc_vend(vc,a,j,Tb);
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < vc_min - cons || vc > vc_max + cons || abs(l-pt) > cons || abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.4')
        end
    
        return;
    end
end

% CASE 4 ------------------ l6 < pt
if(pt > l6)
    vb_upper = vb_max;
    return;
end



error('condition check failed');
end