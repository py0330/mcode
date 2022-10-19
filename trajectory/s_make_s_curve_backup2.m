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
% ------------------ p1 --------------------- 
% 全力减速
% Ta = min(T,s_acc_time(va,     0, a, j))
% Tb = 0
% Tc = T-Ta
% vb = s_acc_vend(va,a,j,Ta)
% v  = vb
% l1 = Ta * (va + vb) / 2.0   
% ------------------ p2 --------------------- 
% 先加速到最大加速度处，再全力减速
% Ta = min(T,2*a/j,s_acc_time(va, max_v, a, j))
% v  = s_acc_vend(va,         a, j,Ta)
% Tb = min(T-Ta, s_acc_time(v, 0, a, j))
% Tc = 0
% vb = s_acc_vend(va + a^2/j,-a,-j,Tb)
% l2 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0
% ------------------ p3 --------------------- 
% 先加速到最大速度处，再全力减速
% Ta = min(T,s_acc_time(va, max_v, a, j))
% v  = s_acc_vend(va, a, j, Ta)
% Tb = min(T-Ta, s_acc_time(v, 0, a, j))
% Tc = 0
% vb = s_acc_vend(v ,-a,-j,Tb)
% l3 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0
% ------------------ p4 --------------------- 
% 先加速到最大速度处，在末端减速到0的情况下，尽可能长的匀速
% Ta = s_acc_time(va, max_v, a, j)
% Tb = min(T-Ta, s_acc_time(max_v, 0, a, j))
% Tc = T-Ta-Tb
% v  = max_v
% vb = 0
% l4 = Ta * (va + v) / 2.0 + Tb * (vb + v) / 2.0 + Tc * max_v
% ------------------ p5 --------------------- 
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
%   v  = va + j*Ta*Ta/4
%   la = j/8*Ta^3 + va*Ta
%   控制距离 p21 : 
%   SUBCASE 1: 可减速到0, 且无匀减速段
%     Tb = 2*sqrt(va/j + Ta*Ta)
%     lb = Tb*v/2
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
%   SUBCASE 2: 可减速到0, 且有匀减速段
%     Tb = v / a + a/j;
%     lb = Tb*v/2
% 
%     带入方程la + lb = pt
%     可得：
%     k4*Ta^4 +k3*Ta^3 + k2*Ta^2 +k1*Ta + k0
%     其中：
%     k4 = j^2/160
%     k3 = j/8
%     k2 = ((j*(va/5 + 5/j))/8 + (j*va)/40)
%     k1 = va
%     k0 = (va*(va/5 + 5/j))/2 - pt
%
%   SUBCASE 3: 不可减速到0, 且无匀减速段
%     Tb = T-Ta
%     lb = Tb*v/2
%
%     带入方程la + lb = pt
%     可得：
%     k2*Ta^2 +k1*Ta + k0
%     其中：
%     k2 = (T*j)/8
%     k1 = va/2
%     k0 = (T*va)/2 - pt
%
%   SUBCASE 4: 不可减速到0, 且有匀减速段
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

% % T1 = 2*a/j
% % T2 = s_acc_time(va,    v2, a, j)
% % T3 = s_acc_time(va, max_v, a, j)
% % T4 = s_acc_time(va,     0, a, j)
% % T5 = s_acc_time(va,    vb, a, j)
% % T6 = s_acc_time(max_v,  0, a, j)
% % T7 = s_acc_time(max_v, vb, a, j)
% % T8 = s_acc_time(v8   ,  0, a, j)
% %
% % v2 = a^2/j
% % v8 = va + a^2/j



% CASE 1: a段减速直到时间结束或速度为0，此时行进的p为最小可能的距离 p1
%   Ta = min(T,T1) , Tb = 0   , Tc = T-Ta
%   p1  = va + 


%
%   在以上情况下时，p > pt，行进距离必定超出pt


% CASE 1: 最大行进长度小于目标长度
%   若 T <  T3 时，Ta = T , Tb = 0   , Tc = 0
%   或 T >= T3 时，Ta = T3, Tb = T-T3, Tc = 0
%
%   在以上情况下时，p < pt，行进距离必定小于pt
%
% CASE 2: 最小行进长度大于目标长度
%   若 T <  T4 时，Ta = T , Tb = 0   , Tc = 0
%   或 T >= T4 时，Ta = T4, Tb = 0   , Tc = T-T4
%
%   在以上情况下时，p > pt，行进距离必定超出pt
%
% CASE 3: 末端行进速度，无法满足vb的要求
%   若 T <  T5 && va > vb时
%
%   在以上情况下时，p > pt，末端速度无法降到vb
%
% CASE 4: vb可以为0
%   若 T >= T1 + T6时，Ta = T1, Tb= T6, Tc = T-T1-T6
%   若 T4 <= T < T1 + T6 时：
%      subcase 1: a&b段均有匀加速：T>=T1+T8
%      subcase 2: a 无匀加速段，b有匀加速段:T1+T2<T<T1+T8
%      subcase 3: a&b段均无匀加速:T<T1+T2
%
% CASE 5: vb可以为0
%   若 T >= T1 + T6时，Ta = T1, Tb= T6, Tc = T-T1-T6
%   若 T4 <= T < T1 + T6 时：
%      subcase 1: a&b段均有匀加速：T>=T1+T8
%      subcase 2: a 无匀加速段，b有匀加速段:T1+T2<T<T1+T8
%      subcase 3: a&b段均无匀加速:T<T1+T2
%
%        % 总体过程应分三种情况：
        % v1 v2 均有匀加速段
        % v1 无匀加速段，v2有匀加速段
        % v1 v2 均无匀加速段 
        
        % CASE 1:
        % v1 v2 处于正好有匀加速段时，
        % 此时应有：v-v1=a^2/j
        % t1 = (v-v1) / a + a/j
        % t2 = (v-v2) / a + a/j
        % T_ = t1 + t2 = (2/a)*v + (2*a)/j - v1/a - v2/a
        %    = (4*a)/j + v1/a - v2/a
        %
        % 若T > T_，则正好有匀加速段
        %
        % CASE 2:
        % v1 无匀加速段，v2有匀加速段
        % v = v2+a^2/j
        % t1 = 2 * sqrt((v-v1) / j )
        % p1 = (v+v1)/2*t1
        % t2 = (v-v2) / a + a/j
        % p2 = (v+v2)/2*t2
        % T_ = t1 + t2
        %    = 2*((a^2/j - v1 + v2)/j)^(1/2) + (2*a)/j
        % 若v1 < v2+a^2/j 且 T > T_，则处于case2
        %
        % CASE 3:
        % 其余情况
%
%
%






% CASE 1：最大速度可以达到 max_v，有能力减速到0，且减速到0时长度满足要求
%   条件：Ta_max + Tb_max < T
%      && pacc + pdec + pmid >= pt
%   此时：vb = 0
%
%          ---                 max_v
%        /     \
% va ---        \
%                \
%                  --- vb      0
%
% CASE 2：最大速度可以达到 max_v，有能力减速到0，但减速到0时长度无法满足要求
%   条件：Ta_max + Tb_max < T
%      && pacc_max + pdec_max + pmid < pt
%   此时：
%     vb = max_v - j*Tb*Tb/4                  |    a/j >  Tb/2
%          max_v - Tb*a + a^2/j               |    a/j <= Tb/2
%     pb = (vb + max_v)/2*Tb
%        = Tb*(max_v - (Tb^2*j)/8)            |    a/j >  Tb/2
%          Tb*(a^2/(2*j) - (Tb*a)/2 + max_v)  |    a/j <= Tb/2
%     pmid = (T - Ta_max - Tb)*max_v
%     p  = pacc_max + pb + pmid
%        = pacc_max - (j*Tb^3)/8 + max_v*(T - Ta)                  |    a/j >  Tb/2
%          (Tb*a^2)/(2*j) - (Tb^2*a)/2 + pacc_max + max_v*(T - Ta) |    a/j <= Tb/2
%
%   在p = pt作为方程时，分别可求得：
%     Tb = sqrt3((pacc_max + max_v*(T - Ta) - pt)*8/j)
%         -B+sqrt(B^2-4C)/2       : B = -a/j C=- (2*pacc_max - 2*pt + 2*max_v*(T - Ta))/a
%               
%
%   其中 p 在 Tb = 2a/j 处的值为：
%     pacc_max - a^3/j^2 + max_v*(T - Ta)
%               
%   因而上述两条件可用 pacc_max - a^3/j^2 + max_v*(T - Ta_max) > pt 代替
%               
%   带入Tb，最终可以求得vb
%          ---                 max_v
%        /     \
% va ---        \
%                \
%                  --- vb      不为0
%
% CASE 3：最大速度无法达到 max_v，此时达到的最高速为 v
%   条件：(Ta_max + Tb_max < T
%        && pacc_max + pdec_max + pmid > pt)
%      or
%        (Ta_max + Tb_max >= T
%        && pacc_max + pdec > pt)                
%
%   此时，判断四个条件：
%     condition 1: Ta >= 2*a/j && Ta >= T-2*a/j  （有匀加速，无匀减速）
%     condition 2: Ta >= 2*a/j && Ta <  T-2*a/j  （有匀加速，有匀减速）
%     condition 3: Ta <  2*a/j && Ta >= T-2*a/j  （无匀加速，无匀减速）
%     condition 4: Ta <  2*a/j && Ta <  T-2*a/j  （无匀加速，有匀减速）
%
%   以上四个条件，所行进的长度依次减少
%   四个条件对应的最短行进长度与最长行进长度所对应的Ta为：
%   [max(2*a/j,T-2*a/j), Tacc_max                    ]
%   [2*a/j             , min(Tacc_max,T-2*a/j)       ]
%   [max(0, T-2*a/j)   , min(Tacc_max, 2*a/j)        ]
%   [0                 , min(Tacc_max, 2*a/j,T-2*a/j)]
%
%   令Ta1 = 2*a/j; Ta2 = T-2*a/j; Ta3 = Tacc_max;
%   对应的最大速度为：v1 v2 v3
%   v1  = s_acc_vend(va,a,j,Ta1)
%   v2  = s_acc_vend(va,a,j,Ta2)
%   v3  = s_acc_vend(va,a,j,Ta3)
%
%   对应的最大末端速度为：
%   vb1 = s_acc_vend(v1,-a,-j,T-Ta1)
%   vb2 = s_acc_vend(v2,-a,-j,T-Ta2)
%   vb3 = s_acc_vend(v3,-a,-j,T-Ta3)
%
%   对应行进的距离为：
%   p1 = (va + v1)/2*Ta1 + (v1 + vb1)/2*(T-Ta1)
%   p2 = (va + v2)/2*Ta2 + (v2 + vb2)/2*(T-Ta2)
%   p3 = (va + v3)/2*Ta3 + (v3 + vb3)/2*(T-Ta3)
%   令计算Ta = 0时行进距离：
%   p0 = (v1 + vb1)/2*T
%
%   原始的四个条件可以描述为：
%   condition 1:
%   [(Ta1 <  Ta2) && (Ta2 < Ta3) && p2 < pt < p3] || 
%   [(Ta1 >= Ta2) && (Ta1 < Ta3) && p1 < pt < p3]
%
%   condition 2:
%   [(Ta2 <  Ta3) && (Ta1 < Ta2) && p1 < pt < p2] || 
%   [(Ta2 >= Ta3) && (Ta1 < Ta3) && p1 < pt < p3]
%
%   condition 3:
%   [(0 <  Ta2) && (Ta1 < Ta2) && p1 < pt < p2] || 
%   [(0 <  Ta2) && (Ta1 < Ta3) && p1 < pt < p3] ||
%   [(0 >= Ta2) && (Ta1 < Ta2) && p1 < pt < p2] || 
%   [(0 <  Ta2) && (Ta1 < Ta3) && p1 < pt < p3] ||
%
%     vb   = v - j*Tb*Tb/4                  |    a/j >  Tb/2
%            v - Tb*a + a^2/j               |    a/j <= Tb/2
%     pacc = (va + v)/2*Ta
%          = Ta*(v - (Ta^2*j)/8)            |    a/j >  Tb/2
%            Ta*(a^2/(2*j) - (Ta*a)/2 + v)  |    a/j <= Tb/2
%     pdec = (vb + v)/2*Tb
%          = Tb*(v - (Tb^2*j)/8)            |    a/j >  Tb/2
%            Tb*(a^2/(2*j) - (Tb*a)/2 + v)  |    a/j <= Tb/2
%     p    = pacc + pdec
%          = pacc - (j*Tb^3)/8 + v*Tb   |    a/j >  Tb/2
%            pacc + (a^2/(2*j) + v)*Tb - (a*Tb^2)/2  |    a/j <= Tb/2
%
%   在p = pt作为方程时，分别可求得：
%     Tb = sqrt3((pacc_max - pt)*8/j)
%         -B+sqrt(B^2-4C)/2       : B = -a/j C=- (2*pacc_max - 2*pt + 2*max_v*(T - Ta))/a
%               
%
%               其中 p 在 Tb = 2a/j 处的值为：
%                  pacc_max - a^3/j^2 + max_v*(T - Ta)
%               
%               因而上述两条件可用 pacc_max - a^3/j^2 + max_v*(T - Ta_max) > pt 代替
%               
%               带入Tb，最终可以求得vb
Ta_max = s_acc_time(va, max_v, a, j);
Tb_max = s_acc_time(0, max_v, a, j);
pacc_max = Ta_max*(va+max_v)/2;
pdec_max = Tb_max*(0+max_v)/2;

% 二分法找vb，尽可能让结束时的速度最小
if(s_curve_length(va, vb_min, max_v, a, j, T) > pt)
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


% 计算满足加减速时间要求的 v 的range
% 可能包含0\1\2个区间:
% CASE 1: 0 个区间
% CASE 2: 1 个区间
%         此时，vmax2 > vmin1
% CASE 3: 2 个区间
%          ---                 vmax1
%        /     \
% va ---        \
%                \
%                  --- vb      
%
% va ---
%        \       
%          ---                 vmin1  
%              \
%                \
%                  --- vb
%          ---                 vmax1
%        /     \
% va ---        \
%                \
%                  --- vb      
%
% 
% va ---
%        \       
%          \                    
%            ---               vmax2
%                \
%                  --- vb
% va ---
%        \       
%         \        --- vb
%          \     /
%            ---               vmin2
%
% vrange 返回 [vmax1, vmin1, vmax2, vmin2]这样的数组，可能为0/2/4维
function vrange = s_curve_vrange(va, vb, a, j, T)
    cons = eps * 10;

    % CHECK CASE 1
    if(s_acc_time(va,vb,a,j) > T)
        vrange = [];
        return;
    end

    v_large = max(va, vb);
    
    % 计算 vmax，vmax 是大于 v_large 的最大可行速度
    v_upper = v_large + a*T;
    v_below = v_large;

    while(abs(v_upper - v_below)>cons)
        v_next = (v_upper + v_below)/2;
        if(s_acc_time(va,v_next,a,j) + s_acc_time(vb,v_next,a,j) > T)
            v_upper  = v_next;
        else
            v_below = v_next;
        end
    end
    vmax = v_below;
    
    % vmin 和 vmax 对称，但需考虑0的影响
    v_small = min(va, vb);
    vmin = max(0, v_small - (vmax - v_large));
    
    % 判断是CASE 2还是3
    vmid = (va+vb)/2;
    if(s_acc_time(va,vmid,a,j) + s_acc_time(vb,vmid,a,j) > T)
        % CASE 3
        v_upper = v_large;
        v_below = vmid;
    
        while(abs(v_upper - v_below)>cons)
            v_next = (v_upper + v_below)/2;
            if(s_acc_time(va,v_next,a,j) + s_acc_time(vb,v_next,a,j) > T)
                v_below = v_next;
            else
                v_upper = v_next;
            end
        end
        vrange = [vmax, v_upper, v_large - v_upper + v_small, vmin];
    else
        % CASE 2
        vrange = [vmax, vmin];
    end

end



