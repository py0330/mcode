function [vb, vc, Ta, Tb, mode] = s_make_s_curve(pa, va, pb, vc_max, a, j, T)
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
%
% 计算满足va vb的曲线
% 需分2种情况考虑：
% ----------------------------------------------------
% CASE A：v > max(va,vb)
%
%           --------          v
%         /          \
% va ----             \
%                      \
%                        ---- vb      
%
% |<- Ta ->|<- Tc ->|<- Tb ->| 
%
% 此时，真实的v可能会无法取到，因而real_v可能会改变. 
% 若 v > max(va,vb), 则real_v 在会取[max(va,vb) v] 中的最大值
% -----------------------------------------------------
% CASE B：v在va,vb之间，例如va > v > vb
%
% va -------
%             \       
%                \          
%                   \        
%                      ------- vb
%  |<- Ta ->|<- Tc ->|<- Tb ->|
% 此时，轨迹仅有1段降速过程，即从va降速到vb，
% 但在va和vb处，会各有一段匀速过程
% 这两段匀速过程持续的时间 按va-v 和 vb-v的比例来决定
%
% 
cons = 100*eps;

va = max(va,0);

Z1 = a^2/j;
Z2 = T^2*j;

pt = pb - pa;
T_va_to_max_v = s_acc_time(va, vc_max, a, j);
T_0_to_va     = s_acc_time(va, 0, a, j);

% ------------------ l1 --------------------- %
% 全力减速
Ta = min(T, T_0_to_va);
vb = max(s_acc_vend(va,-a,-j,Ta), 0);
l  = Ta * (va + vb) / 2.0;
if(pt < l - 1e-10)
    % CASE 1
    va
    a
    j
    T
    s_curve_length(va, 0, 0, a, j, T)
    pt
    
    [Tmax,Tmin]=s_s_curve_Tmax_Tmin( ...
            pa, ...
            va, ...
            pb, ...
            3, ...
            vc_max, ...
            a, ...
            j)

    error('s_make_s_curve failed: too SMALL distance');
end

% ------------------ l2 --------------------- %
% vb为0，无需加速即可完成
% 此时 0 < v < va, 因此Ta段匀速va，Tb段匀速0，Tc段减速到0 
Tc = T_0_to_va;
Ta = T - Tc;
l  = Ta * va + Tc * va / 2;
if(pt < l)
    Tc = T_0_to_va;
    Ta = (pt - Tc * va / 2)/va;
    Tb = T - Ta - Tc;
    vc  = Ta / (Tb + Ta) * va;
    vb = 0;
    mode = 1;
    % error %
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 || ...
            Ta > T + cons || Ta < -cons ||...
            Tb > T + cons || Tb < -cons ||...
            Tc > T + cons || Tc < -cons)
        Ta
        Tb
        vc
        vb
        s_curve_length(va, 0, vc, a, j, T)
        pt

        error('l2 error');

    end


    return;
end

% ------------------ l3 --------------------- %
% vb为0，a段达不到最大加速度，b段达不到最大加速度
%
% 此时需要满足3个条件：
% A. Ta + Tb <= T
% B. v - va  <= a^2/j && v <= a^2/j && v <= max_v && v >= va
% C. pt <= l
%
% 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
% Ta = 2 * sqrt((v-va)/j)
% Tb = 2 * sqrt(v/j)
% 带入条件 A，以下为推导对 v 的表达式：
% eq = (Ta+Tb)^2
% expand(eq)
% >> 8*(v/j - va/j)^(1/2)*(v/j)^(1/2) + (8*v)/j - (4*va)/j
% left  = 8*(v/j - va/j)^(1/2)*(v/j)^(1/2)
% right = T^2 + (4*va)/j - (8*v)/j
%
% 原不等式等效于 left^2 - right^2 <= 0 && right >= 0
%
% expand(left^2 - right^2)
% collect(left^2 - right^2,v)
% >> ((16*(T^2 + (4*va)/j))/j - (64*va)/j^2)*v - (T^2 + (4*va)/j)^2
% solve(left^2 - right^2,v)
% >> (T^4*j^2 + 8*T^2*j*va + 16*va^2)/(16*T^2*j)
% 对于不等式 left^2 - right^2 <= 0，应有：
% v <= (T^4*j^2 + 8*T^2*j*va + 16*va^2)/(16*T^2*j)
% && v <= (j*T^2)/8 + va/2
% 考虑其他条件应有：
% v_upper = min(max_v, a^2/j 
%               ,(T^4*j^2 + 8*T^2*j*va + 16*va^2)/(16*T^2*j)
%               ,(j*T^2)/8 + va/2)
% v_below = va
v_upper = min([vc_max, Z1, (Z2^2 + 8*Z2*va + 16*va^2)/(16*Z2),Z2/8 + va/2]);
v_below = va;
l = -1;
if(v_upper >= v_below)
    vc = v_upper;
    Ta = s_acc_time(va,vc,a,j);
    Tb = s_acc_time(0, vc,a,j);
    Tc = T - Ta - Tb;
    l = (va + vc)/2*Ta + vc/2*Tb + vc*Tc;
end
if(pt < l)
    % syms va j Ta pt T f(Ta)
    % v  = va + j*Ta*Ta/4
    % la = (v+va)/2*Ta
    % Tb = sqrt(4*va/j + Ta*Ta)
    % lb = Tb*v/2
    % lc = (T - Ta - Tb)*v
    %
    % l = la+lb+lc
    %   = la + T*v - Ta*v -Tb*v/2
    Ta = newton_raphson_binary_search(@(Ta)(...
        T*va - (Ta^3*j)/8 ...
        - (va/2 + Ta^2*j/8)*(Ta^2 + 4*va/j)^(1/2) ...
        + (T*Ta^2*j)/4 - pt)...
        ,0,min([s_acc_time(va,vc,a,j),T/2,T-s_acc_time(0,va,a,j)])...
        ,10*eps);
    vc  = va + j*Ta*Ta/4;
    vb = 0;
    Tb = sqrt(4*va/j + Ta*Ta);
    mode = 0;
    % error %
    if(abs(s_curve_length(va, 0, va + j*Ta*Ta/4, a, j, T) - pt) > 1e-10)
        Ta = newton_raphson_binary_search(@(Ta)(...
            (va + j*Ta*Ta/8)*Ta...
            + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
            + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4) - pt)...
            ,0,min([s_acc_time(va,vc,a,j),T/2,T-s_acc_time(0,va,a,j)])...
            ,10*eps);
        error('error');
    end
    
    
    return;
end

% ------------------ l4 --------------------- %
% vb为0，a段达不到最大加速度，b段可达到最大加速度
%
% 此时需要满足3个条件：
% A. Ta + Tb <= T
% B. v - va  <= a^2/j && v >= a^2/j && v <= max_v && v >= va
% C. pt <= l
%
% 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
% Ta = 2 * sqrt((v-va)/j)
% Tb = v/a + a/j
% 带入条件 A, 有：
% Ta^2 - (T-Tb)^2 <= 0 && Tb <= T
% 展开有：
% -v^2 + k1*v + k0 <=0 && v  <= T*a-a^2/j 
% 其中：
% k1 = (2*T*a + 2*a^2/j)
% k0 = -a^2*((T - a/j)^2 + (4*va)/j)
%
% 可以看出，其极值点为 a^2/j + T*a，然而若v 大于此值，必有 Tb > T
% 故而取其较小的根 x1
% 因此有 v <= x1
% 进一步的，根据条件 B，可以得到v的取值范围：
%
% v_upper = min(x1, va + a^2/j, max_v, T*a-a^2/j)
% v_below = a^2/j
k1 = 2*(T*a + Z1);
k0 = -a^2*((4*va)/j + T^2 + Z1/j - (2*T*a)/j);
x1 = (k1-sqrt(k1*k1+4*k0))/2;% k0为-1

v_upper = min([x1, va + Z1, vc_max, T*a-Z1]);
v_below = max(Z1, va);

if(v_upper < v_below)
    l=-1;
else
    vc = v_upper;
    Ta = s_acc_time(va,vc,a,j);
    Tb = s_acc_time(0, vc,a,j);
    Tc = T - Ta - Tb;
    l = (va + vc)/2*Ta + vc/2*Tb + vc*Tc;
end
if(pt < l)
    % clear
    % syms va j Ta pt T a f(Ta)
    % v  = va + j*Ta*Ta/4
    % la = (v+va)/2*Ta
    % Tb = (v + a^2/j)/a
    % lb = Tb*v/2
    % lc = (T - Ta - Tb)*v
    %
    % l = la+lb+lc
    %
    % collect(l-pt,Ta)
    %
    % k4 = -j^2/32
    % k3 = -j/8
    % k2 = (j*(a^2/j + va))/8 - (j*va)/8 - (j*(a^2/j - T + va))/4
    % k0 = (va*(a^2/j + va))/2 - va*(a^2/j - T + va) - pt
    %
    % 考虑4次方程过于复杂，这里用迭代方法
    k4 = -j^2/(32*a);
    k3 = -j/8;
    k2 = (T*j)/4 - a/8 - (j*va)/(4*a);
    k0 = T*va - pt - va^2/(2*a) - (a*va)/(2*j);
    
    % 求其上下界
    % syms f(Ta) g(Ta)
    % f(Ta) = k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k0
    % g(Ta) = 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta
    % solve(g,Ta)
    %
    % 可得 f 的三个极值点：
    % r0 = 0
    % r1 = (-3*a + (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j);
    % r2 = (-3*a - (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j);
    % 
    % 分析极值点分布，有：
    % r1 = (-3*a + (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j)
    %    = (-3*a + (a^2 + 16*j*(T*a-va))^(1/2))/(2*j)
    %
    % 由于b段必定可达加速度a，因此有 T >= Tb = v/a+a/j >= va/a+a/j
    % 因此有：
    % r1 =  (-3*a + (a^2 + 16*j*(T*a-va))^(1/2))/(2*j)
    %    >= (-3*a + (17*a^2)^(1/2))/(2*j)
    %    >  0
    % 而显然 r2 < 0
    % 因此其极值点分布为：
    % r2 < r0 = 0 < r1
    %
    % 下求 f(0)：
    % f(0) =  T*va - pt - va^2/(2*a) - (a*va)/(2*j)
    %      =  -pt + va*(Ta + Tb + Tc - (va/a + a/j)/2)
    %      =  -pt + va * Ta + va * (Tb-(va/a + a/j)/2) + va * Tc
    %      <= -pt + (va+v)/2*Ta + va*(v/a-va/(2*a)+a/(2*j)) + v*Tc
    %      =  -pt + la + lc + 1/(2*a)*(2*va*v - va*va + va*a*a/j)
    %      =  -pt + la + lc + 1/(2*a)*(v*v-(v*v-2*v*va +va*va) + va*a*a/j)
    %      =  -pt + la + lc + 1/(2*a)*(v*v + va*a*a/j)
    %      <= -pt + la + lc + 1/(2*a)*(v*v + v*a*a/j)
    %      =  -pt + la + lc + v*(v/a + a/j)/2
    %      =  -pt + la + lc + lb
    %      =  0
    %
    % 因此可用 newton_raphson 方法搜索 [0,r1] 内的取值
    Ta_below = 0;
    Ta_upper = (-3*a + (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j);
    Ta = newton_raphson_binary_search(@(x)(k4*x^4 + k3*x^3 + k2*x^2 + k0),Ta_below,Ta_upper,cons);

    % return
    vc  = va + j*Ta*Ta/4;
    vb = 0;
    Tb = (vc + Z1)/a;
    mode = 0;
    % error
    if(abs(s_curve_length(va, 0, vc, a, j, T) - pt) > 1e-10)
        s_curve_length(va, 0, vc, a, j, T);
        error('error in l4');
    end
    
    
    return;
end

% ------------------ l5 --------------------- %
% vb为0，a段达不到最大加速度，b段可达到最大加速度
%
% 此时需要满足3个条件：
% A. Ta + Tb <= T
% B. v - va  >= a^2/j && v >= a^2/j && v <= max_v && v >= va
% C. pt <= l
%
% 首先根据条件 A 和 B 确定 v 的取值范围，然后再计算 l
% Ta = (v-va)/a + a/j
% Tb = v/a + a/j
% 带入条件 A, 有：
% v <= - a^2/j + (T*a)/2 + va/2
%
% v_upper = min(max_v, - a^2/j + (T*a)/2 + va/2)
% v_below = va + a^2/j
v_upper = min(vc_max, - Z1 + (T*a)/2 + va/2);
v_below = va + Z1;
l = -1;
if(v_below <= v_upper)
    vc = v_upper;
    Ta = (vc-va)/a + a/j;
    Tb = vc/a + a/j;
    Tc = T - Ta - Tb;
    l = Ta*(vc+va)/2 + Tb*vc/2 + Tc * vc;
end
if(pt < l)
    % clear
    % syms va j Ta pt T a f(Ta)
    % v  = va + Ta*a - a^2/j
    % la = (v+va)/2*Ta
    % Tb = (v + a^2/j)/a
    % lb = Tb*v/2
    % lc = (T - Ta - Tb)*v
    %
    % l  = la+lb+lc-pt
    %
    % collect(l,Ta)
    %
    % k2 = -a
    % k1 = (a*(T - va/a) + a^2/j)
    % k0 = (T - va/a)*(va - a^2/j) + (va*(va - a^2/j))/(2*a)- pt
    
    k2 = -a;
    k1 = (a*T - va + Z1);
    k0 = (va - Z1)*(T - va/a/2)- pt;

    Ta = (-k1 + sqrt(k1*k1 - 4*k0*k2))/(2*k2);

    % return
    vc  = va + Ta*a - Z1;
    vb = 0;
    Tb = (vc + Z1)/a;
    mode = 0;
    % error
    if(abs(s_curve_length(va, 0, vc, a, j, T) - pt) > 1e-10)
        s_curve_length(va, 0, vc, a, j, T)
        error('error in l5');
    end
        
    return;
end

% ------------------ l6 --------------------- % 
% vb不为0，达不到max_v，a段达不到最大加速度，b段可达到最大加速度
%
% 此时需要满足3个条件：
% A. Ta + Tb =  T
% B. v - va  <= a^2/j && v - vb >= a^2/j && v <= max_v && v >= va
% C. pt <= l
%
% 条件 A   可得 Ta  = T - Tb
%            => Ta <= T - 2*a/j   (因为 Tb >= 2*a/j)
%
% 条件 B.1 可得 Ta <= 2*a/j
% 条件 B.2 可得 va + j*Ta*Ta/4 >= a^2/j
%            => Ta >= 2*sqrt(max(0, a^2/j - va)/j)
% 条件 B.3 可得 Ta <= T_va_to_max_v
% 条件 B.4 可得 Ta >= 0  (包含在B.2中)
Ta_upper = min([T_va_to_max_v, 2*a/j, T - 2*a/j]);
Ta_below = 2*sqrt(max(0, Z1 - va)/j);
if(Ta_upper >= Ta_below)
    Ta = Ta_upper;
    vc  = va + j*Ta*Ta/4;
    la = Ta*(vc+va)/2;
    Tb = T-Ta;
    vb = vc - Tb*a + Z1;
    lb = Tb*(vc+vb)/2;
    l = la + lb;
end
if(pt < l)
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

    vc  = va + j*Ta*Ta/4;
    Tb = T - Ta;
    vb = s_acc_vend(vc,-a,-j,Tb);
    mode = 0;
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
            || Tb<-cons || Ta < -cons)
        error('wrong Tb in l6')
    end
    return;
end

% ------------------ l8 -------------------- 
% vb不为0，达不到max_v，a段可达到最大加速度，b段可达到最大加速度
Ta_upper = min([T_va_to_max_v, T - 2*a/j]);
Ta_below = 2*a/j;
l=-1;
if(Ta_upper >= Ta_below)
    Ta = Ta_upper;
    vc  = va + Ta*a - Z1;
    la = (va + vc)*Ta/2;
    Tb = T-Ta;
    vb = s_acc_vend(vc,-a,-j,Tb);
    lb = Tb*(vc+vb)/2;
    l  = la + lb;
end
if(pt < l + cons)
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
    vb = s_acc_vend(vc,-a,-j,Tb);
    mode = 0;
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
            ||Tb<-cons)
        error('wrong Tb in l8')
    end
    return;
end

% ------------------ l10 -------------------- 
% vb不为0，可达max_v，b段可达到最大加速度
Ta = min(s_acc_time(va, vc_max, a, j),T);
vc  = vc_max;
la = (va + vc_max)*Ta/2;
Tb = min([T - Ta, s_acc_time(0, vc_max, a, j), 2*a/j]);
if(Tb < 2*a/j)
    l = -1;
else
    vb = s_acc_vend(vc,-a,-j,Tb);
    lb = Tb*(vc+vb)/2;
    Tc = T - Ta - Tb;
    lc = Tc * vc_max;
    l  = la + lb + lc;
end
if(pt < l + cons)
    B = -a/j;
    C = -(2*la - 2*pt + 2*vc_max*(T - T_va_to_max_v))/a;

    Tb = max((-B+sqrt(B^2-4*C))/2,0);
    Ta = s_acc_time(va, vc_max, a, j);
    vb = vc_max - Tb*a + Z1;
    vc  = vc_max;
    mode = 0;
    % error %
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
            ||Tb<-cons || Ta < -cons)
        error('wrong Tb in l10')
    end
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 以下对应可能产生的最大l


% ------------------ l7 --------------------- 
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
Ta_upper = min([T, T_va_to_max_v, 2*a/j]);
Ta_below = max(0, T - 2*a/j);
l = -1; % 这里为必要条件，因为有可能 vb = v
if(Ta_upper - Ta_below >= 0)
    Ta = Ta_upper;
    vc  = va + j*Ta*Ta/4;
    Tb = T-Ta;
    vb = max(vc - j*Tb*Tb/4, 0);
    l  = Ta * (va + vc) / 2.0 + Tb * (vb + vc) / 2.0;
end
% 以下条件1可能存在边界点的误判，例如只有减速段为0时
% 因此若一定无法达到匀加速状态且一定无法加速到max_v，则必然进入该条件
if(pt < l || T < min(2*a/j, T_va_to_max_v))
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
    k1 = -3*Z2/8;
    k0 = pt - (T*(2*va - Z2/4))/2;

    % 选根
    % 其极值为 (k1)/(2*k2) = (3*T)/2
    % 因此需选其较小的根
    Ta = (-k1-sqrt(k1*k1-4*k0*k2))/2/k2;
    
    vc  = va + j*Ta*Ta/4;
    Tb = T - Ta;
    vb = s_acc_vend(vc,-a,-j,Tb);
    mode = 0;
    % error %
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
            ||Tb < -cons)
        error('wrong Tb in l7')
    end
    return;
end

% ------------------ l9 -------------------- 
% vb不为0，达不到max_v，a段可达到最大加速度，b段不可达到最大加速度
Ta_upper = min([T, T_va_to_max_v]);
Ta_below = 2*a/j;
if(Ta_upper >= Ta_below)
    Ta = Ta_upper;
    vc  = s_acc_vend(va,a,j,Ta);
    la = (va + vc)*Ta/2;
    Tb = T-Ta;
    vb = s_acc_vend(vc,-a,-j,Tb);
    lb = Tb*(vc+vb)/2;
    l  = la + lb;
end
if(pt < l || T < T_va_to_max_v)
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
    mode = 0;
    % error %
    if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
            || Tb<-cons || Ta < -cons)
        Ta = -1;
        r = cubic_equation_solve3(k3,k2,k1,k0);
        for i=1:3
            if(isreal(r(i)) && r(i) > -1e-10)
                Ta = abs(r(i));
                break;
            end
        end
        error('wrong Tb in l9')
    end
    return;
end


% ------------------ l11 -------------------- 
% vb不为0，可达max_v，b段不可达到最大加速度
Ta = T_va_to_max_v;
la = (va + vc_max)*Ta/2;
Tb = 0;
lb = 0;
Tc = T - Ta - Tb;
lc = Tc * vc_max;
l  = la + lb + lc;
% 此处最终先check error
if(pt >= l + cons)
    error('s_make_s_curve failed: too LARGE distance');
    return;
end
Tb = max((la + vc_max*(T - T_va_to_max_v) - pt)*8/j,0)^(1/3);
vb = vc_max - j*Tb*Tb/4;
vc  = vc_max;
mode = 0;
if(abs(s_curve_length(va, vb, vc, a, j, T) - pt) > 1e-10 ...
        || Tb < -cons || Ta < -cons || Ta + Tb > T + cons ...
        || vb < -cons ...
        || abs((va + vc)*Ta/2 + (vb + vc)*Tb/2 + vc*(T - Ta - Tb) - pt) > 1e-10)
    Ta
    Tb
    Tc

    (va + vc)*Ta/2 + (vb + vc)*Tb/2 + vc*Tc
    pt
    s_curve_length(va, vb, vc, a, j, T)

    error('wrong Tb in l11')
end





end





