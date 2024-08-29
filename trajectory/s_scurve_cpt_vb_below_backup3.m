function [vb_below] = s_scurve_cpt_vb_below(pa, va, pb, vc_max, vb_max, a, j, T)
% 在给定过程中起始速度va, 最大速度v，加速度a，跃度j，时间长度T的情况下
% 自动计算末端所可能达到的【最大的】vb
%
% t      : current time
% pa     : init pos
% va     : init vel
% pb     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period
%
cons = 100*eps;
pt = pb - pa;

vb_below = s_scurve_cpt_vb1(va, pt, a, j, T, vb_max, vc_max);
if(vb_below >= 0)
    return;
end

vb_below = s_scurve_cpt_vb2(va, pt, a, j, T, vb_max, vc_max);
if(vb_below >= 0)
    return;
end

vb_below = s_scurve_cpt_vb3(va, pt, a, j, T, vb_max, vc_max);
if(vb_below >= 0)
    return;
end

end


% CASE 1: vb = 0
function [vb_below] = s_scurve_cpt_vb1(va, pt, a, j, T, vb_max, vc_max)
    % 判断是否属于当前情况，即可以达到vb_upper
    % CASE 1.1 : vc = 0, Tc > 0
    % CASE 1.2 : v1 <= vc <= v2, Tc = 0   (v1 = min(va, vb_max), v2 = max(va, vb_max))
    % CASE 1.3 : vc <= v2, Tc = 0   (v1 = max(va, vb_max), v2 = min(va, vb_max))
    %      1.3.1 : va -> vc 无匀速段，vc -> vb 无匀速段
    %      1.3.2 : va -> vc 无匀速段，vc -> vb 有匀速段
    %      1.3.3 : va -> vc 有匀速段，vc -> vb 无匀速段
    %      1.3.4 : va -> vc 有匀速段，vc -> vb 有匀速段
    
    % v1  是 min(va,vb)
    % v2  是 max(va,vb)
    % ve1 是 v1 - a*a/j
    % ve2 是 v2 - a*a/j
    % 
    % 正常的大小关系应该为：
    % ve1 < ve2 < v1 < v2
    % 或 
    % ve1 < v1 < ve2 < v2
    %
    % 因此，以上case的条件为：
    % v1         < vc < v2          : CASE 1.2
    % max(ve2,0) < vc < v1          : CASE 1.3.1
    % max(ve1,0) < vc < min(ve2,v1) : CASE 1.3.2 OR 1.3.3
    % 0          < vc < ve1         : CASE 1.3.4
    %              vc = 0           : CASE 1.1 
    % 因此对于 T 来说，应该有：
    % T_va_to_vb < T                             : CASE 1.2
    % T_va_to_vb < T < T_va_to_ve2 + T_vb_to_ve2 : CASE 1.3.1
    % ve2 < v1 ? T_va_to_vb : T_va_to_ve2 + T_vb_to_ve2
    %       < T < T_va_to_ve1 + T_vb_to_ve1      : CASE 1.3.2 OR 1.3.3
    % T_va_to_ve1 + T_vb_to_ve1 < T < T_va_to_0 + T_vb_to_0     
    %                                            : CASE 1.3.4
    % T_va_to_0 + T_vb_to_0 < T                  : CASE 1.1
    
    Z1 = a^2/j;
    Z2 = T^2*j;

    v1 = min(va,0);
    v2 = max(va,0);
    ve1 = v1 + a*a/j;
    ve2 = v2 + a*a/j;
    
    T_va_to_vb    = s_acc_time(va, 0, a, j);
    T_va_to_vcmax = s_acc_time(va, vc_max, a, j);
    T_va_to_ve1   = s_acc_time(va, ve1, a, j);
    T_va_to_ve2   = s_acc_time(va, ve2, a, j);
    T_vb_to_vcmax = s_acc_time(0, vc_max, a, j);
    T_vb_to_ve1   = s_acc_time(0, ve1, a, j);
    T_vb_to_ve2   = s_acc_time(0, ve2, a, j);

    %------------------------ CASE 1.1 -------------------------%
    % CASE 1.1 : vc = 0, Tc > 0
    % 
    % 若要减速到0，需要满足:
    % A. T >= T_va_to_vcmax + T_vb_to_vcmax
    % B. l = T_va_to_vcmax*(va+vcmax)/2 + T_vb_to_vcmax*(vc_max)/2 <= pt
    if(T_va_to_vcmax + T_vb_to_vcmax <=  T)
        l = T_va_to_vcmax * (va+vcmax) / 2 + T_vb_to_vcmax * vc_max / 2;
        if(l >= pt)
            vb_below = 0;

            % debug check %
            if(T_va_to_vcmax * (va+vcmax)/2 + T_vb_to_vcmax * vc_max/2 < pt)
                error('wrong vb_below in CASE 1.1')
            end

            return;
        end
    end
    
    %------------------------ CASE 1.2 -------------------------%
    % CASE 1.2 : v1 <= vc <= v2, Tc = 0
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T
    % B. l = T_va_to_vb*(va + vb_max)/2 + (T - T_va_to_vb) * vb_max <= pt
    if(T >= T_va_to_vb)
        l = T_va_to_vb * (v1+v2)/2 + (T-T_va_to_vb)*v2;
        if(l >= pt)
            vb_below = 0;
            
            % debug check %
            if(T_va_to_vb*(va + vb_below)/2 + (T - T_va_to_vb) * va < pt)
                error('wrong vb_upper in CASE 1.2')
            end

            return;
        end      
    end
    
    %------------------------ CASE 1.3.1 -----------------------%
    % CASE 1.3.1 : va -> vc 无匀速段，vc -> vb 无匀速段
    % 
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T <= T_va_to_ve1 + T_vb_to_ve1
    % B. l = Ta*(va + vc)/2 + Tb * (0 + vc)/2 <= pt
    %
    % 其中 vc 的范围是 [v2, min(ve1, vc_max)]
    %
    % 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
    % Ta = 2 * sqrt((vc-va)/j)
    % Tb = 2 * sqrt((vc-vb)/j)
    % 带入条件 Ta + Tb == 0，可得：
    % vc_ans = solve(Ta + Tb == T, vc)
    % >> (T^4*j^2 + 8*T^2*j*va + 8*T^2*j*vb + 16*va^2 - 32*va*vb + 16*vb^2)/(16*T^2*j)
    %
    vc_below = v2;
    vc_upper = min(ve1,vc_max);
    T_below  = T_va_to_vb;
    T_upper  = T_va_to_ve1 + T_vb_to_ve1;
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)
        
        vb = 0;
        vc = (T^4*j^2 + 8*T^2*j*va + 8*T^2*j*vb + 16*va^2 - 32*va*vb + 16*vb^2)/(16*T^2*j);
        Ta = s_acc_time(va,vc,a,j);
        Tb = s_acc_time(vb, vc,a,j);
        l = Ta*(va+vc)/2 + Tb*(0+vc)/2;
        if(l >= pt)
            vb_below = 0;

            % debug check %
            if(vc < 0 || vc > vc_max)
                error('wrong vb_below in CASE 1.3.1')
            end

            return;
        end
    end

    %------------------------ CASE 1.3.2 & 1.3.3 -----------------------%
    % CASE 1.3.2 : va -> vc 无匀速段，vc -> vb 有匀速段
    % CASE 1.3.3 : va -> vc 有匀速段，vc -> vb 无匀速段
    % 
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T <= T_va_to_vcmax + T_vb_to_vcmax
    % B. l = Ta*(va + vc)/2 + Tb * (0 + vc)/2 <= pt
    %
    % 其中 vc 的范围是 [0, min(v2 - a^2/j, 0)]
    %
    % 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
    % T1 = 2 * sqrt((vc-v1)/j)
    % T2 = (vc-v2)/a + a/j
    % 带入条件 A，以下为推导对 vc 的表达式：
    % solve(T1 + T2 == T, vc)
    % >> (j*vb_max - a^2 + 2*a*j*((va - vb_max + T*a)/j)^(1/2) - T*a*j)/j
    %    -(a^2 - j*vb_max + 2*a*j*((va - vb_max + T*a)/j)^(1/2) + T*a*j)/j
    %
    % ve2 < v1 ? T_va_to_vb : T_va_to_ve2 + T_vb_to_ve2
    %       < T < T_va_to_ve1 + T_vb_to_ve1      : CASE 1.3.2 OR 1.3.3
    vc_below = max(ve1, v2);
    vc_upper = min(ve2, vc_max);
    T_below  = s_acc_time(v1,vc_below,a,j) + s_acc_time(v2,vc_below,a,j);
    T_upper  = s_acc_time(v1,vc_upper,a,j) + s_acc_time(v2,vc_upper,a,j);
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)
        
        % 选根 tbd
%         vc = (a^2 + T*a*j + 2*a*j*(-(va - T*a)/j)^(1/2))/j;
        vc = (a^2 + T*a*j - 2*a*j*(-(va - T*a)/j)^(1/2))/j;
        T1 = s_acc_time(v1,vc,a,j);
        T2 = s_acc_time(v2, vc,a,j);
        l = T1*(v1+vc)/2 + T2*(v2+vc)/2;
        if(l >= pt)
            vb_below = 0;

            % debug check %
            if(vc < 0 || vc > vc_max)
                error('wrong vb_upper in CASE 1.3.1')
            end

            return;
        end
    end

    %------------------------ CASE 1.3.4 -----------------------%
    % CASE 1.3.4 : va -> vc 有匀速段，vc -> vb 有匀速段
    % 
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T <= T_va_to_0 + T_0_to_vb
    % B. l = Ta*(va + vc)/2 + Tb * (vb_max + vc)/2 <= pt
    %
    % 其中 vc 的范围是 [0, min(v1 - a^2/j, 0)]
    %
    % 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
    % T1 = (v1-vc)/a + a/j
    % T2 = (v2-vc)/a + a/j
    % 带入条件 A，以下为推导对 vc 的表达式：
    % solve(T1 + T2 == T, vc)
    % >> (j*vb_max - a^2 + 2*a*j*((va - vb_max + T*a)/j)^(1/2) - T*a*j)/j
    %    -(a^2 - j*vb_max + 2*a*j*((va - vb_max + T*a)/j)^(1/2) + T*a*j)/j
    %
    vc_below = ve2;
    vc_upper = vc_max;
    T_below  = T_va_to_ve2 + T_vb_to_ve2;
    T_upper  = T_va_to_vcmax + T_vb_to_vcmax;
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)
        
        vc = (- 2*a^2 + T*j*a + j*v1 + j*v2)/(2*j);
        T1 = s_acc_time(v1, vc,a,j);
        T2 = s_acc_time(v2, vc,a,j);
        l = T1*(v1+vc)/2 + T2*(v2+vc)/2;
        if(l >= pt)
            vb_below = 0;

            % debug check %
            if(vc < 0 || vc > vc_max)
                error('wrong vb_upper in CASE 1.3.1')
            end

            return;
        end
    end
    
    vb_below = -1;
    return;
end

% CASE 2: 无匀速段，vb = [0, vb_max]
function [vb_below] = s_scurve_cpt_vb2(va, pt, a, j, T, vb_max, vc_max)
    % CASE 2.1: va -> vc 无匀加速 vc -> vb 无匀加速
    % CASE 2.2: va -> vc 无匀加速 vc -> vb 有匀加速
    % CASE 2.3: va -> vc 有匀加速 vc -> vb 有匀加速
    % CASE 2.4: va -> vc 有匀加速 vc -> vb 无匀加速

    cons = 100*eps;
    Z1 = a^2/j;
    Z2 = T^2*j;

    %------------------------ CASE 2.1 -----------------------%
    % CASE 2.1: va -> vc 无匀加速 vc -> vb 无匀加速
    % vb不为0，达不到max_v，a段达不到最大加速度，b段达不到最大加速度
    %
    % 此时需要满足3个条件：
    % A. Ta + Tb =  T
    % B. v - va  <= a^2/j && v - vb <= a^2/j && v <= max_v && v >= va
    % C. pt <= l
    %
    % 条件 A   可得 Ta =  T-Tb
    %            => Ta >= T - 2*a/j   (因为Tb < 2*a/j)
    %
    % 条件 B.1 可得 Ta <= 2*a/j
    % 条件 B.2 无法得到有效等式，因为vb可以为为任意值
    % 条件 B.3 可得 Ta <= T_va_to_max_v
    % 条件 B.4 可得 Ta >= 0
    T_va_to_vcmax = s_acc_time(va, vc_max, a, j);
    Ta_upper = min([T, 2*a/j, T_va_to_vcmax]);
    Ta_below = max(0, T - 2*a/j);
    l = -1;
    if(Ta_below <= Ta_upper)
        Ta = Ta_upper;
        vc = va + j*Ta*Ta/4;
        Tb = T-Ta;
        vb = max(vc - j*Tb*Tb/4, 0);
        l_upper = Ta * (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;

        Ta = Ta_below;
        vc = va + j*Ta*Ta/4;
        Tb = T-Ta;
        vb = max(vc - j*Tb*Tb/4, 0);
        l_below = Ta * (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
    end
    if(pt <= l_upper && pt >= l_below)
        % v  = va + j*Ta*Ta/4
        % la = j/8*Ta^3 + va*Ta
        % Tb = T-Ta
        % lb = Tb*(v+vb)/2
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
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.1')
        end

        return;
    end
    
    %------------------------ CASE 2.2 -----------------------%
    % vb不为0，达不到max_v，a段达不到最大加速度，b段可达到最大加速度
    %
    % 此时需要满足3个条件：
    % A. Ta + Tb =  T
    % B. vc - va <= a^2/j && vc - vb >= a^2/j && v <= max_v && v >= va &&
    % vb > 0
    % C. pt <= l
    %
    % 条件 A   可得 Ta  = T - Tb
    %            => Ta <= T - 2*a/j   (因为 Tb >= 2*a/j)
    %
    % 条件 B.1 可得 Ta <= 2*a/j
    % 条件 B.2 可得 va + j*Ta*Ta/4 - vb >= a^2/j
    %            => j*Ta*Ta/4 >= a^2/j - va + vb
    %            => j*Ta*Ta/4 >= a^2/j - va
    %            => Ta >= 2*sqrt(max(0, a^2/j - va)/j)
    % 条件 B.3 可得 Ta <= T_va_to_max_v
    % 条件 B.4 可得 Ta >= 0  (包含在B.2中)
    % 条件 B.5 可得 vb = va + j*Ta*Ta/4 + a^2/j - a*(T-Ta) >= 0
    %            => j*Ta*Ta/4 - a*Ta + a*T + va - a^2/j >= 0
    %            => Ta > 
    l = -1;
    A = j/4;
    B = a;
    C = -a*T + va + a^2/j;
    if(B*B-4*A*C >=0)
        Ta_below2 = (-B + sqrt(B*B-4*A*C))/(2*A);
    else
        Ta_below2 = 0
    end

    Ta_upper = min([T_va_to_vcmax, 2*a/j, T - 2*a/j]);
    Ta_below = max(2*sqrt(max(0, a*a/j-va)/j), Ta_below2);
    if(Ta_upper >= Ta_below)
        Ta = Ta_upper;
        vc  = va + j*Ta*Ta/4;
        la = Ta*(vc+va)/2;
        Tb = T-Ta;
        vb = vc - Tb*a + a*a/j;
        lb = Tb*(vc+vb)/2;
        l_upper = la + lb;

        Ta = Ta_below;
        vc  = va + j*Ta*Ta/4;
        la = Ta*(vc+va)/2;
        Tb = T-Ta;
        vb = vc - Tb*a + a*a/j;
        lb = Tb*(vc+vb)/2;
        l_below = la + lb;


    end
    if(pt <= l_upper && pt >= l_below)
        %   clear
        %   syms va j Ta a T
        %   v  = va + j*Ta*Ta/4
        %   la = j/8*Ta^3 + va*Ta
        %   Tb = T-Ta
        %   vb = v - Tb*a + a^2/j;
        %   lb = Tb*(v+vb)/2
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
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.2')
        end

        return;
    end

    %------------------------ CASE 2.3 -----------------------%
    % CASE 2.3: va -> vc 有匀加速 vc -> vb 有匀加速
    %
    Ta_upper = min([T_va_to_vcmax, T - 2*a/j]);
    Ta_below = 2*a/j;
    l=-1;
    if(Ta_upper >= Ta_below)
        Ta = Ta_upper;
        vc  = va + Ta*a - Z1;
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,-a,-j,Tb);
        lb = Tb*(vc+vb)/2;
        l_upper  = la + lb;

        Ta = Ta_below;
        vc  = va + Ta*a - Z1;
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,-a,-j,Tb);
        lb = Tb*(vc+vb)/2;
        l_below  = la + lb;
    end
    if(pt <= l_upper && pt >= l_below)
        % syms va T a j pt Ta
        % v  = va + Ta*a - a^2/j
        % la = Ta*(va + v)/2
        % Tb = T-Ta
        % vb = v - Tb*a + a^2/j;
        % lb = Tb*(v+vb)/2
        %
        % 根据 la + lb = pt，有：
        % collect(la+lb-pt,Ta)
        % - a*Ta^2 + 2*T*a*Ta - pt - (T*(T*a - 2*va + a^2/j))/2
        %
        % 可得方程系数
        k2 = -a;
        k1 = 2*T*a;
        k0 = - pt - (T*(Z1 + T*a - 2*va))/2;

        % 选根
        % 其极值为：
        % r = k1 / (2*k2) = T
        % 应有 Ta < T
        % 故而选其较小的根
        Ta = (-k1+sqrt(k1*k1-4*k0*k2))/2/k2;

        vc  = va + Ta*a - Z1;
        Tb = T - Ta;
        vb_below = s_acc_vend(vc,-a,-j,Tb);
        vb = vb_below;

        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_below in CASE 2.3')
        end

        return;
    end

    %------------------------ CASE 2.4 -----------------------%
    % CASE 2.4: va -> vc 有匀加速 vc -> vb 无匀加速
    %
    % 此时需要满足3个条件：
    % A. Ta + Tb =  T
    % B. va - vc <= a^2/j && vb - vc <= a^2/j && vc >= 0 && vc <= va
    % C. pt >= l
    %
    % 条件 A   可得 Ta =  T - Tb
    %            => Ta >= T - 2*a/j   (因为Tb < 2*a/j)
    %
    % 条件 B.1 可得 Ta <= 2*a/j
    % 条件 B.2 无法得到有效等式，因为vb可以为为任意值
    % 条件 B.3 可得 Ta <= T_va_to_0
    % 条件 B.4 可得 Ta >= 0
    Ta_upper = min([T, T_va_to_vcmax]);
    Ta_below = 2*a/j;
    l = -1; % 这里为必要条件，因为有可能 vb = v
    if(Ta_upper >= Ta_below)
        Ta = Ta_upper;
        vc  = s_acc_vend(va,a,j,Ta);
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,-a,-j,Tb);
        lb = Tb*(vc+vb)/2;
        l_upper  = la + lb;

        Ta = Ta_below;
        vc  = s_acc_vend(va,a,j,Ta);
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,-a,-j,Tb);
        lb = Tb*(vc+vb)/2;
        l_below = la + lb;
    end
    % 以下条件1可能存在边界点的误判，例如只有减速段为0时
    % 因此若一定无法达到匀加速状态且一定无法加速到max_v，则必然进入该条件
    if(l >= 0 && pt >= l)
        % syms va T a j pt Ta
        % v  = va + Ta*a - a^2/j
        % la = Ta*(va + v)/2
        % Tb = T-Ta
        % vb = v - j*Tb*Tb/4;
        % lb = Tb*(v+vb)/2
        %
        % 根据 la + lb = pt，有：
        % collect(la+lb-pt,Ta)
        % (j*Ta^3)/8 + (- a/2 - (3*T*j)/8)*Ta^2 + ((T^2*j)/8 + a^2/(2*j) + (T*(2*a + (T*j)/2))/2)*Ta - pt - (T*((T^2*j)/4 - 2*va + (2*a^2)/j))/2
        %
        % 可得方程系数
        k3 = j/8;
        k2 = - a/2 - (3*T*j)/8;
        k1 = Z2*3/8 + Z1/2 + T*a;
        k0 = - pt - (T*(Z2/4 + 2*Z1 - 2*va))/2;

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

        Ta = newton_raphson_binary_search(@(x)(k3*x*x*x+k2*x*x+k1*x+k0),2*a/j,T,10*eps);

        vc  = va + Ta*a - Z1;
        Tb = max(T - Ta, 0);
        vb = s_acc_vend(vc,-a,-j,Tb);
        vb_below = vb;
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons || abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.4')
        end

        return;
    end

    vb_below = -1;

end

% CASE 3: 有匀速段，vc = vc_max, vb > 0
function [vb_below] = s_scurve_cpt_vb3(va, pt, a, j, T, vb_max, vc_max)
    % CASE 3.1: vc -> vb 有匀加速
    % CASE 3.2: vc -> vb 无匀加速
    
    cons = 100*eps;

    % vc = 0;
    Ta = s_acc_time(va, vc_max, a, j);
    la = Ta*(va+vc_max)/2;

    T_vbmax_to_vcmax = s_acc_time(vb_max, vc_max, a, j);
    T_0_to_vcmax = s_acc_time(0, vc_max, a, j);

    % CASE 3.1:
    Tb_upper = min(T-Ta, T_0_to_vcmax);
    Tb_below = max(2*a/j, T_vbmax_to_vcmax);
    l = -1;
    if(Tb_upper >= Tb_below)
        Tb = Tb_below;
        vb = s_acc_vend(vc_max,-a,-j,Tb);
        lb = Ta*(vb+vc_max)/2;
        l  = la + lb + (T-Ta-Tb)*vc_max;
    end
    if(l >= 0 && l >= pt)
        B = -a/j;
        C = -(2*la - 2*pt + 2*vc_max*(T - Ta))/a;

        Tb = max((-B+sqrt(B^2-4*C))/2,0);
        Ta = s_acc_time(va, vc_max, a, j);
        vb = vc_max - Tb*a + Z1;

        vb_below = vb;

        % debug check %
        if(abs(Tb*vb/2 + Ta*va/2 + (T-Ta-Tb)*vc_max -pt) > cons)
            error('wrong vb_upper in CASE 3.1')
        end

        return;
    end


    % CASE 3.2:
    Tb_upper = min([T-Ta, 2*a/j, T_vbmax_to_vcmax]);
    Tb_below = 0;
    l = -1;
    if(Tb_upper >= Tb_below)
        Tb = Tb_below;
        vb = s_acc_vend(vc_max,-a,-j,Tb);
        lb = Tb*(vb+vc_max)/2;
        l  = la + lb + (T-Ta-Tb)*vc_max;
    end
    if(l >= 0 && l >= pt)
        Tb = max((la + vc_max*(T - Ta) - pt)*8/j,0)^(1/3);
        vb = vc_max - j*Tb*Tb/4;
        
        vb_below = vb;

        % debug check %
        if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
            error('wrong vb_upper in CASE 3.2')
        end

        return;
    end
    

    vb_below = -1;
    return;
end