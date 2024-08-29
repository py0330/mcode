function [vb_upper] = s_scurve_cpt_vb_upper(pa, va, pb, vc_max, vb_max, a, j, T)
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

vb_upper = s_scurve_cpt_vb1(va, pt, a, j, T, vb_max, vc_max);
if(vb_upper >= 0)
    return;
end

vb_upper = s_scurve_cpt_vb2(va, pt, a, j, T, vb_max, vc_max);
if(vb_upper >= 0)
    return;
end

vb_upper = s_scurve_cpt_vb3(va, pt, a, j, T, vb_max, vc_max);
if(vb_upper >= 0)
    return;
end

end


% CASE 1: vb = vb_max
function [vb_upper] = s_scurve_cpt_vb1(va, pt, a, j, T, vb_max, vc_max)
    % 判断是否属于当前情况，即可以达到vb_upper
    % CASE 1.1 : vc = 0, Tc > 0
    % CASE 1.2 : v1 <= vc <= v2, Tc = 0   (v1 = max(va, vb_max), v2 = min(va, vb_max))
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

    v1 = min(va,vb_max);
    v2 = max(va,vb_max);
    ve1 = v1 - a*a/j;
    ve2 = v2 - a*a/j;
    
    T_va_to_vb  = s_acc_time(va, vb_max, a, j);
    T_va_to_0   = s_acc_time(va, 0, a, j);
    T_va_to_ve1 = s_acc_time(va, ve1, a, j);
    T_va_to_ve2 = s_acc_time(va, ve2, a, j);
    T_vb_to_0   = s_acc_time(vb_max, 0, a, j);
    T_vb_to_ve1 = s_acc_time(vb_max, ve1, a, j);
    T_vb_to_ve2 = s_acc_time(vb_max, ve2, a, j);

    %------------------------ CASE 1.1 -------------------------%
    % CASE 1.1 : vc = 0, Tc > 0
    % 
    % 若要减速到0，需要满足:
    % A. T >= T_va_to_0 + T_0_to_vb
    % B. l = T_va_to_0*va/2 + T_0_to_vb*vb_max/2 <= pt
    if(T_va_to_0 + T_vb_to_0 <=  T)
        l = T_va_to_0 * va / 2 + T_vb_to_0 * vb_max / 2;
        if(l <= pt)
            vb_upper = vb_max;

            % debug check %
            if(T_va_to_0*va/2 + T_vb_to_0*vb_upper/2 > pt)
                error('wrong vb_upper in CASE 1.1')
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
        l = T_va_to_vb * (v1+v2)/2 + (T-T_va_to_vb)*v1; % 此时 vc 恰好等于 v1
        if(l <= pt)
            vb_upper = vb_max;
            
            % debug check %
            if(T_va_to_vb*(va + vb_upper)/2 + (T - T_va_to_vb) * min(va,vb_upper) > pt)
                error('wrong vb_upper in CASE 1.2')
            end

            return;
        end      
    end
    
    %------------------------ CASE 1.3.1 -----------------------%
    % CASE 1.3.1 : va -> vc 无匀速段，vc -> vb 无匀速段
    % 
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T <= T_va_to_0 + T_0_to_vb
    % B. l = Ta*(va + vc)/2 + Tb * (vb_max + vc)/2 <= pt
    %
    % 其中 vc 的范围是 [max(v2 - a^2/j, 0), v1]
    %
    % 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
    % Ta = 2 * sqrt((va-vc)/j)
    % Tb = 2 * sqrt((vb-vc)/j)
    % 带入条件 Ta + Tb == 0，可得：
    % vc_ans = solve(Ta + Tb == T, vc)
    % >> (- T^4*j^2 + 8*T^2*j*va + 8*T^2*j*vb - 16*va^2 + 32*va*vb - 16*vb^2)/(16*T^2*j)
    %
    vc_below = max(ve2,0);
    vc_upper = v1;
    T_below  = T_va_to_vb;
    T_upper  = T_va_to_ve2 + T_vb_to_ve2;
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)
        
        vb = vb_max;
        vc = (- T^4*j^2 + 8*T^2*j*va + 8*T^2*j*vb - 16*va^2 + 32*va*vb - 16*vb^2)/(16*T^2*j);
        
        Ta = s_acc_time(va,vc,a,j);
        Tb = s_acc_time(vb_max, vc,a,j);
        l = Ta*(va+vc)/2 + Tb*(vb_max+vc)/2;
        if(l <= pt)
            vb_upper = vb_max;

            % debug check %
            if(vc < 0 || vc > vc_max)
                error('wrong vb_upper in CASE 1.3.1')
            end

            return;
        end
    end

    %------------------------ CASE 1.3.2 & 1.3.3 -----------------------%
    % CASE 1.3.2 : va -> vc 无匀速段，vc -> vb 有匀速段
    % CASE 1.3.3 : va -> vc 有匀速段，vc -> vb 无匀速段
    % 
    % 此时需要满足以下条件：
    % A. T_va_to_vb <= T <= T_va_to_0 + T_0_to_vb
    % B. l = Ta*(va + vc)/2 + Tb * (vb_max + vc)/2 <= pt
    %
    % 其中 vc 的范围是 [0, min(v2 - a^2/j, 0)]
    %
    % 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
    % T1 = 2 * sqrt((v1-vc)/j)
    % T2 = (v2-vc)/a + a/j
    % 带入条件 A，以下为推导对 vc 的表达式：
    % solve(T1 + T2 == T, vc)
    % >> (j*vb_max - a^2 + 2*a*j*((va - vb_max + T*a)/j)^(1/2) - T*a*j)/j
    %    -(a^2 - j*vb_max + 2*a*j*((va - vb_max + T*a)/j)^(1/2) + T*a*j)/j
    %
    % ve2 < v1 ? T_va_to_vb : T_va_to_ve2 + T_vb_to_ve2
    %       < T < T_va_to_ve1 + T_vb_to_ve1      : CASE 1.3.2 OR 1.3.3
    vc_below = max(ve1, 0);
    vc_upper = min(ve2, v1);
    T_below  = s_acc_time(v1,vc_upper,a,j) + s_acc_time(v2,vc_upper,a,j);
    T_upper  = s_acc_time(v1,vc_below,a,j) + s_acc_time(v2,vc_below,a,j);
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)

        vc = -(a*j - j*v2 + a^2 - 2*a*j*((a + v1 - v2)/j)^(1/2))/j;
        T1 = s_acc_time(v1,vc,a,j);
        T2 = s_acc_time(v2, vc,a,j);
        l = T1*(v1+vc)/2 + T2*(v2+vc)/2;
        if(l <= pt)
            vb_upper = vb_max;

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
    vc_below = 0;
    vc_upper = ve1;
    T_below  = T_va_to_ve1 + T_vb_to_ve1;
    T_upper  = T_va_to_0 + T_vb_to_0;
    if(vc_upper >= vc_below ...
        && T >= T_below ...
        && T <= T_upper)
        
        vc = (2*a^2 - T*j*a + j*v1 + j*v2)/(2*j);
        T1 = s_acc_time(v1, vc,a,j);
        T2 = s_acc_time(v2, vc,a,j);
        l = T1*(v1+vc)/2 + T2*(v2+vc)/2;
        if(l <= pt)
            vb_upper = vb_max;

            % debug check %
            if(vc < 0 || vc > vc_max)
                error('wrong vb_upper in CASE 1.3.1')
            end

            return;
        end
    end
    
    vb_upper = -1;
    return;
end

% CASE 2: 无匀速段，vb = [0, vb_max]
function [vb_upper] = s_scurve_cpt_vb2(va, pt, a, j, T, vb_max, vc_max)
    % CASE 2.1: va -> vc 无匀加速 vc -> vb 无匀加速
    % CASE 2.2: va -> vc 无匀加速 vc -> vb 有匀加速
    % CASE 2.3: va -> vc 有匀加速 vc -> vb 有匀加速
    % CASE 2.4: va -> vc 有匀加速 vc -> vb 无匀加速

    cons = 100*eps;
    
    %------------------------ CASE 2.1 -----------------------%
    % CASE 2.1: va -> vc 无匀加速 vc -> vb 无匀加速
    % 此时需要满足以下条件：
    % A. Ta + Tb = T
    % B. va - vc <= a^2/j && vb - vc <= a^2/j && vc <= va && vc >= 0
    % C. l = Ta*(va+vc)/2 + Tb(vb+vc)/2 <= pt
    %
    % 条件 A   可得 Ta =  T-Tb
    %            => Ta >= T - 2*a/j   (因为Tb < 2*a/j)
    %
    % 条件 B.1 可得 Ta <= 2*a/j
    % 条件 B.2 可得 Ta <= T_va_to_0
    % 条件 B.4 可得 Ta >= 0
    T_va_to_0 = s_acc_time(va, 0, a, j);
    Ta_upper = min([T, 2*a/j, T_va_to_0]);
    Ta_below = max(0, T - 2*a/j);
    l = -1;
    if(Ta_below <= Ta_upper)
        Ta = Ta_upper;
        vc = va - j*Ta*Ta/4;
        Tb = T-Ta;
        vb = max(vc + j*Tb*Tb/4, 0);
        l  = Ta * (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
    end
    if(l >= 0 && pt >= l)
        % vc = va - j*Ta*Ta/4
        % la = -j/8*Ta^3 + va*Ta
        % Tb = T-Ta
        % vb = vc + j*Tb*Tb/4
        % lb = Tb*(v+vb)/2
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

        % debug check %
        vb = vb_upper;
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.1')
        end

        return;
    end
    
    %------------------------ CASE 2.2 -----------------------%
    % CASE 2.2: va -> vc 无匀加速 vc -> vb 有匀加速
    %
    % 此时需要满足3个条件：
    % A. Ta + Tb =  T
    % B. va - vc  <= a^2/j && vb - vc >= a^2/j && vc >= 0 && vc <= va
    % C. pt >= l
    %
    % 条件 A   可得 Ta  = T - Tb
    %            => Ta <= T - 2*a/j   (因为 Tb >= 2*a/j)
    %
    % 条件 B.1 可得 Ta <= 2*a/j
    % 条件 B.2 可得 vb - va + j*Ta*Ta/4 >= a^2/j
    %            => j*Ta*Ta/4 >= va - vb + a^2/j
    %            => j*Ta*Ta/4 >= va - vb_max + a^2/j
    %            => Ta >= 2*sqrt(max(0, va - a^2/j)/j)
    % 条件 B.3 可得 Ta <= T_va_to_0
    % 条件 B.4 可得 Ta >= 0  (包含在B.2中)
    l = -1;
    Ta_upper = min([T_va_to_0, 2*a/j, T - 2*a/j]);
    Ta_below = 2*sqrt(max(0, va - vb_max + a^2/j)/j);
%     Ta_below = 0;
    if(Ta_upper >= Ta_below)
        Ta = Ta_upper;
        vc = va - j*Ta*Ta/4;
        la = Ta*(vc+va)/2;
        Tb = T-Ta;
        vb = vc + Tb*a - a*a/j;
        lb = Tb*(vc+vb)/2;
        l  = la + lb;
    end
    if(l >= 0 && pt >= l)
        %   clear
        %   syms va j Ta a T
        %   vc = va - j*Ta*Ta/4
        %   la = j/8*Ta^3 + va*Ta
        %   Tb = T-Ta
        %   vb = vc - Tb*a + a^2/j;
        %   lb = Tb*(v+vb)/2
        %   l  = la + lb
        %
        %   【result】:
        %   带入方程la + lb = pt
        %   可得：
        %   k3*Ta^3 + k2*Ta^2 + k1*Ta + k0
        %   其中：
        %   k3 = j/8
        %   k2 = -(T*j)/4 + a/2
        %   k1 = a^2/(2*j) - T*a
        %   k0 = (T*(-a^2/j + T*a + 2*va))/2 - pt
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
        k2 = -(T*j)/4 + a/2;
        k1 = a^2/(2*j) - T*a;
        k0 = (T*(-a^2/j + T*a + 2*va))/2 - pt;
        
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
        vb_upper = s_acc_vend(vc,a,j,Tb);
    
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.2')
        end

        return;
    end

    %------------------------ CASE 2.3 -----------------------%
    % CASE 2.3: va -> vc 有匀加速 vc -> vb 有匀加速
    %
    Ta_upper = min([T_va_to_0, T - 2*a/j]);
    Ta_below = 2*a/j;
    l=-1;
    if(Ta_upper >= Ta_below)
        Ta = Ta_upper;
        vc  = va - Ta*a + a*a/j;
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,a,j,Tb);
        lb = Tb*(vc+vb)/2;
        l  = la + lb;
    end
    if(l >= 0 && pt >= l)
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
        k2 = a;
        k1 = -2*T*a;
        k0 = - pt + (T*(a*a/j + T*a + 2*va))/2;
    
        % 选根
        % 其极值为：
        % r = k1 / (2*k2) = T
        % 应有 Ta < T
        % 故而选其较小的根
        Ta = (-k1-sqrt(k1*k1-4*k0*k2))/2/k2;
        
        vc  = va - Ta*a + a*a/j;
        Tb = T - Ta;
        vb = s_acc_vend(vc,a,j,Tb);
        vb_upper = vb;

        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons|| abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.3')
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
    Ta_upper = min([T, T_va_to_0]);
    Ta_below = 2*a/j;
    l = -1; % 这里为必要条件，因为有可能 vb = v
    if(Ta_upper - Ta_below >= 0)
        Ta = Ta_upper;
        vc  = s_acc_vend(va,-a,-j,Ta);
        la = (va + vc)*Ta/2;
        Tb = T-Ta;
        vb = s_acc_vend(vc,a,j,Tb);
        lb = Tb*(vc+vb)/2;
        l  = la + lb;
    end
    % 以下条件1可能存在边界点的误判，例如只有减速段为0时
    % 因此若一定无法达到匀加速状态且一定无法加速到max_v，则必然进入该条件
    if(l >= 0 && pt >= l)
        % syms va j Ta a T pt
        % vc = va - Ta*a + a*a/j
        % la = Ta*(va+vc)/2
        % Tb = T-Ta
        % vb = vc + j*Tb*Tb/4
        % lb = Tb*(vc+vb)/2;
        %
        % 根据 la + lb = pt，有：
        % collect(la+lb-pt,Ta)
        %
        % 可得方程系数
        k3 = -j/8;
        k2 = a/2 + (3*T*j)/8;
        k1 = -T*T*j*3/8 - a*a/j/2 - T*a;
        k0 = - pt + (T*(T*T*j/4 + 2*a*a/j + 2*va))/2;
    
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
    
        vc  = va - Ta*a + a*a/j;
        Tb = max(T - Ta, 0);
        vb = s_acc_vend(vc,a,j,Tb);
        vb_upper = s_acc_vend(vc,a,j,Tb);
        % debug check %
        l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
        if(vc < 0 || vc > vc_max || abs(l-pt) > cons || abs(Ta+Tb-T) > cons)
            error('wrong vb_upper in CASE 2.4')
        end

        return;
    end

    vb_upper = -1;

end

% CASE 3: 有匀速段，vc = 0, vb > 0
function [vb_upper] = s_scurve_cpt_vb3(va, pt, a, j, T, vb_max, vc_max)
    % CASE 3.1: vc -> vb 有匀加速
    % CASE 3.2: vc -> vb 无匀加速
    
    cons = 100*eps;

    % vc = 0;
    Ta = s_acc_time(va, 0, a, j);
    la = Ta*va/2;

    T_vb_to_0 = s_acc_time(vb_max, 0, a, j);

    % CASE 3.1:
    Tb_upper = min(T-Ta, T_vb_to_0);
    Tb_below = 2*a/j;
    l = -1;
    if(Tb_upper >= Tb_below)
        Tb = Tb_below;
        vb = s_acc_vend(0,a,j,Tb);
        lb = Ta*vb/2;
        l  = la + lb;
    end
    if(l >= 0 && l <= pt)
        % vb = Tb*a - a*a/j
        % lb = Tb*vb/2 
        %    = pt - la
        %
        % collect(Tb*vb/2 - pt + la, Tb)
        %
        % gives:
        % k2 * Tb^2 + k1 * Tb + k0 == 0
        %
        % where:
        % k2 = a/2
        % k1 = -a*a/j/2
        % k0 = la - pt
        %
        
        k2 = a/2;
        k1 = -a*a/j/2;
        k0 = la - pt;

        % 上式极值为 -k1/k2/2 = a/j/2, 因此取其大根

        Tb = (-k1 + sqrt(k1*k1-4*k2*k0))/(2*k2);
        vb_upper = s_acc_vend(0,a,j,Tb);

        % debug check %
        if(abs(Tb*vb_upper/2 + Ta*va/2 -pt) > cons)
            error('wrong vb_upper in CASE 3.1')
        end

        return;
    end


    % CASE 3.2:
    Tb_upper = min([T-Ta, 2*a/j, T_vb_to_0]);
    Tb_below = 0;
    l = -1;
    if(Tb_upper >= Tb_below)
        Tb = Tb_below;
        vb = s_acc_vend(0,a,j,Tb);
        lb = Ta*vb/2;
        l  = la + lb;
    end
    if(l >= 0 && l <= pt)
        % vb = j*Tb*Tb/4
        % lb = Tb*vb/2 
        %    = pt - la
        %
        % collect(Tb*vb/2 - pt + la, Tb)
        %
        % gives:
        % k3 * Tb^3 + k0 == 0
        %
        % where:
        % k3 = j/8
        % k0 = la - pt
        %
        k3 = j/8;
        k0 = la - pt;

        Tb = (-k0/k3)^(1/3);
        vb_upper = s_acc_vend(0,a,j,Tb);

        % debug check %
        if(abs(Tb*vb_upper/2 + Ta*va/2 -pt) > cons)
            error('wrong vb_upper in CASE 3.2')
        end

        return;
    end
    

    vb_upper = 0;
    return;
end