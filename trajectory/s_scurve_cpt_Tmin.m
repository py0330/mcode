function Tmin = s_scurve_cpt_Tmin(pa, va, pb, vb_max, vc_max, a, j)
% 计算当前点位所需的最大最小时间
%
% pa     : init pos
% va     : init vel
% pb     : end  pos
% max_vb : max  end vel
% v      : max  vel  during period
% a      : max  acc  during period
% j      : max  jerk during period
% T      : period
%
% Tmax：开始时尽可能快的减速，若减速到0，则为inf，否则以到达pb的时间为准
% Tmin：开始时尽可能快的加速，直到速度最大，之后保持最大速度到终点

cons = eps * 10000;
pt = pb - pa;
Z1 = a^2/j;

T_va_to_vb = s_acc_time(va,vb_max,a,j);
l_va_to_vb = T_va_to_vb*(va + vb_max) /2;
if(va > vb_max && l_va_to_vb > pt + max(pt,1)*1e-10)
    Tmin=-1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 计算Tmin：%%%%%%

% ------------------ l1 --------------------
% 加速不到max_vb, 无法达到最大加速度a
% va < vb < max_vb, vb - va < a^2 / j
vb = min(va + Z1, vb_max);
l = min(s_acc_time(va,vb,a,j) * (va + vb)/2, l_va_to_vb);
if(va < vb_max && pt < l)
%     %%%%%%%%%%%%%%%%% METHOD1 %%%%%%%%%%%%%%%%% 
%     % clear
%     % syms va vb a j pt
%     % T = 2 * sqrt((vb-va) / j )
%     % l = T * (va+vb)/2
%     % expand((T*(va + vb)/2)^2 *j - pt^2*j)
%     % coeffs(l^2*j - pt^2*j, vb)
%     %
%     % k3 = 1;
%     % k2 = va;
%     % k1 = -va^2;
%     % k0 = -j*pt^2 - va^3;
%     %
%     % f(vb) = k3*vb^3 + k2*vb^2 + k1*vb + k0 
%     % 
%     % 对f求导，可知
%     % df    = 3*k3*vb^2 + 2*k2*vb + k1
%     %
%     % 其极值为 -va 与 1/3va
%     % 将此2值带入，可以看出
%     % f(-va)  = -j*pt^2
%     % f(va/3) = - j*pt^2 - (32*va^3)/27
%     % 因此该三次方程在pt不为0时只有一个根，为0时有3个根，此时选择最大的
%     k3 = 1;
%     k2 = va;
%     k1 = -va^2;
%     k0 = -j*pt^2 - va^3;
%     
%     r = cubic_equation_solve(k3,k2,k1,k0);
%     % 选根 %
%     % 【注意】 pt为0时有3个根，两个为 -va，因此这里选择最大的
%     if(isreal(r(3)))
%         vb = r(3);
%     else
%         vb = r(1);
%     end
    
    %%%%%%%%%%%%%%%%% METHOD2 %%%%%%%%%%%%%%%%% 
    % newton raphson %
    vb = newton_raphson_binary_search(@(x)(sqrt((x-va)/j)*(va+x) - pt)...
        , va, vb_max, 10*eps);
    Tmin = s_acc_time(va,vb,a,j);
    return;
end

% ------------------ l2 --------------------
% 加速不到max_vb, 可以达到最大加速度a
% va < vb < max_vb, vb - va >= a^2 / j
l = l_va_to_vb;
if(va < vb_max && pt < l)
    % clear
    % syms va vb a j pt;
    % T = (vb-va) / a + a/j;
    % l = T * (va+vb)/2;
    % coeffs(l,vb)
    % 
    % k2 = 1/(2*a);
    % k1 = a/(2*j);
    % k0 = (va*(a/j - va/a))/2 - pt;
    % 
    % 其极值为 vb = -k1/2/k2 = -a^2/(2*j)
    % 因此取其右侧大值

    k2 = 1/(2*a);
    k1 = a/(2*j);
    k0 = (va*(a/j - va/a))/2 - pt;

    % 选根 %
    vb = (-k1 + sqrt(k1*k1-4*k2*k0))/(k2*2);
    Tmin = s_acc_time(va,vb,a,j);
    return;
end

% ------------------ l3 --------------------
% 可以加速到某个v，无匀速段，之后减速到max_vb，加减速过程，加速度均不超过a
% v2 = max(va,vb), v1 = min(va,vb), v2 < v < min(v1 + a^2/j, max_v)
vb = vb_max;
v2 = max(va,vb);
v1 = min(va,vb);
v_upper  = min(v1 + Z1, vc_max);
v_below  = v2;
l = -1;
if(v_upper >= v_below)
    l = s_acc_time(v2, v_upper, a, j) * (v_upper + v2)/2 +...
        s_acc_time(v1, v_upper, a, j) * (v_upper + v1)/2;
end
if(pt < l)
    % 以下为方程求解
    % syms v2 v1 a j pt v;
    % T1 = 2 * sqrt((v-v2) / j )
    % T2 = 2 * sqrt((v-v1) / j )
    % l = T1 * (v2 + v) / 2 + T2 * (v1 + v) / 2
    % 
    % 有：
    % diff(l,v) == T1/2 + T2/2 + (v + v2)/(j*T1) + (v + v1)/(j*T2)
    % newton-raphson method & binary search
    v = newton_raphson_binary_search(@(x)(...
        sqrt((x-v2)/j)*(v2+x) + sqrt((x-v1)/j)*(v1+x) - pt)...
        , v_below, v_upper, 10*eps);

    T1 = 2 * sqrt((v-v2) / j );
    T2 = 2 * sqrt((v-v1) / j );
    Tmin = T1 + T2;
    return;
end

% ------------------ l4 --------------------
% 可以加速到某个v，无匀速段，之后减速到max_vb，加减速过程，有一段可以达到加速度a
% v2 = max(va,vb), v1 = min(va,vb), v2 < v < min(v1 + a^2/j, max_v)
v2 = max(va,vb);
v1 = min(va,vb);
v_upper  = min(v2 + Z1, vc_max);
v_below  = max(v1 + Z1, v2);
l = -1;
if(v_upper >= v_below)
    l = s_acc_time(v2, v_upper, a, j) * (v_upper + v2)/2 +...
        s_acc_time(v1, v_upper, a, j) * (v_upper + v1)/2;
end
if(pt < l)
    % 以下为方程求解
    % syms v2 v1 a j pt v;
    % T1 = 2 * sqrt((v-v2) / j )
    % T2 = (v-v1) / a + a/j
    % l1 = T1 * (v2 + v) / 2
    % l2 = T2 * (v1 + v) / 2
    % l  = l1 + l2
    % 
    % 可以化成1元4次方程：
    % eq = expand(l1^2 - (pt-l2)^2)
    % 
    % 受限于难以求解1元4次方程，因此还是使用牛顿法 
    %
    % eq = l - pt
    % >> eq  = sqrt((x-v2)/j)*(v2+x) + ((x-v1)/a+a/j)*(v1+x)/2 - pt)
    % deq = diff(eq,v)
    % >> deq = T1/2 + T2/2 + (v + v1)/(2*a) + (v+v2)/(j*T1)
    % 
    % vc = newton_raphson_binary_search(@(x)(...
    %         sqrt((x-v2)/j)*(v2+x) + ((x-v1)/a+a/j)*(v1+x)/2 - pt)...
    %         , v_below, v_upper, 10*eps);
    %
    % 因为 vc 在接近 v2 附近时，上式数值性能不佳，因此改用T1做未知数
    T1 = newton_raphson_binary_search(@(T1)(...
        T1*(v2+v2+j*T1*T1/4)/2 + ((v2+j*T1*T1/4-v1)/a+a/j)*(v1+v2+j*T1*T1/4)/2 - pt)...
        , 0, 2*a/j, 10*eps);
    vc = v2+j*T1*T1/4;
    T2 = (vc-v1) / a + a/j;
    Tmin = T1 + T2;
    return;
end

% ------------------ l5 --------------------
% 可以加速到某个v，无匀速段，之后减速到max_vb，加减速过程，两段都可以达到加速度a
% v2 = max(va,vb), v1 = min(va,vb), v2 + a^2/j < v < max_v
v2 = max(va,vb);
v1 = min(va,vb);
v_upper  = vc_max;
v_below  = v2 + Z1;
if(v_upper < v_below)
    l = -1;
else
    l = s_acc_time(v2, v_upper, a, j) * (v_upper + v2)/2 +...
        s_acc_time(v1, v_upper, a, j) * (v_upper + v1)/2;
end
if(pt < l)
    % 以下为方程求解
    % syms v2 v1 a j pt v;
    % T1 = (v-v2) / a + a/j
    % T2 = (v-v1) / a + a/j
    % l1 = T1 * (v2 + v) / 2
    % l2 = T2 * (v1 + v) / 2
    % l  = l1 + l2
    % 
    % 可以化成1元2次方程：
    % f(v) = l-pt
    % f(v) == 0
    % 
    % 其系数为：
    % k2 = 1/a
    % k1 = a/j
    % k0 = (v2*(a/j - v2/a))/2 - pt + (v1*(a/j - v1/a))/2
    %
    % 分析 f 的极值，有：
    % df = diff(f,v)
    % solve(df,v)
    %
    % ans = -a^2/(2*j)
    %
    % 因此f的根位于上述结果两侧，应取较大值
    k2 = 1/a;
    k1 = a/j;
    k0 = (v2*(a/j - v2/a))/2 - pt + (v1*(a/j - v1/a))/2;
    
    v  = (-k1+sqrt(k1*k1 - 4*k0*k2))/(2*k2);

    T1 = (v-v2) / a + a/j;
    T2 = (v-v1) / a + a/j;
    Tmin = T1 + T2;
    return;
end

% ------------------ l6 --------------------
% 可以加速到max_v，匀速运行一段时间，之后减速到max_vb
% v2 = max(va,vb), v1 = min(va,vb), v = max_v
v2 = max(va,vb);
v1 = min(va,vb);
v  = vc_max;
T1 = s_acc_time(v2,vc_max,a,j);
T2 = s_acc_time(v1,vc_max,a,j);
T3 = (pt - T1 * (v + v2)/2 - T2 * (v + v1)/2)/vc_max;
Tmin = T1 + T2 + T3;
return;

end
