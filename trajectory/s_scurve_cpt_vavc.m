function [va,vc,Ta,Tb,mode] = s_scurve_cpt_vavc(pa, pb, vb, va_upper, va_below, vc_max, a, j, T)
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

Ta = 0;
Tb = 0;
mode = 0;


cons = 1000*eps;
pt = pb - pa;

% 根据速度可达修正 va_range
va_upper = min(va_upper, s_acc_vend(vb, a, j, T));
va_below = max(va_below, s_acc_vend(vb, -a, -j, T));

% CASE 1: vc > max(vb,va_range)
% va = va_upper
% vc > max(va_upper, vb)
T_va_upper_to_vb  = s_acc_time(va_upper, vb, a, j);
l_upper = max(va_upper, vb) * (T - T_va_upper_to_vb) + T_va_upper_to_vb*(vb+va_upper)/2;
if(pt > l_upper)
    % 如果用 vc 做未知数，用 newton-ranphson 方法，在vc 很接近v2的时候数值性能很差
    % 例如：
    va = va_upper;
    vc_upper = s_cpt_vc_upper_by_va_vb_T(va,vb,T,a,j);

    vc = newton_raphson_binary_search(@(vc)(...
        s_acc_time(va,vc,a,j)*(va+vc)/2  ...
        + s_acc_time(vb,vc,a,j)*(vb+vc)/2  ...
        + max(T-s_acc_time(va,vc,a,j)-s_acc_time(vb,vc,a,j), 0)*vc ...
        - pt) ...
        ,max(va,vb),min(vc_max, vc_upper)...
        ,10*eps);
    Ta = s_acc_time(va,vc,a,j);
    Tb = s_acc_time(vb,vc,a,j);
    mode = 0;
    l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0 + (T-Ta-Tb)*vc;
    
%     va = va_upper;
%     v1 = min(va,vb);
%     v2 = max(va,vb);
%     
%     % 事实上 newton-ranphson 内的方程应该为：
%     % vc = s_acc_vend(v2,a,j,T2);
%     % T1 = s_acc_time(v1,vc,a,j);
%     %     
%     % T1*(v1+vc)/2 + T2*(v2+vc)/2 + (T-T2-T1)*vc - pt
%     
%     % 但matlab似乎只支持单行函数
%     T2_upper = min(newton_raphson_binary_search(@(T2)(...
%         s_acc_vend(v2,a,j,T2) - s_acc_vend(v1,a,j,T-T2)) ...
%         ,0,T...
%         ,eps), s_acc_time(v2,vc_max,a,j));
%     
%     T2 = newton_raphson_binary_search(@(T2)(...
%         s_acc_time(v1,s_acc_vend(v2,a,j,T2),a,j)*(v1+s_acc_vend(v2,a,j,T2))/2  ...
%         + T2*(v2+s_acc_vend(v2,a,j,T2))/2  ...
%         + (T-T2-s_acc_time(v1,s_acc_vend(v2,a,j,T2),a,j))*s_acc_vend(v2,a,j,T2) ...  % 本行不同
%         - pt) ...
%         ,0,T2_upper...
%         ,10*eps);
% 
% 
%     vc = s_acc_vend(v2,a,j,T2);
%     T1 = s_acc_time(v1,vc,a,j);
%     if(va < vb)
%         Ta = T1;
%         Tb = T2;
%     else
%         Tb = T1;
%         Ta = T2;
%     end


    % debug check %
    l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0 + (T-Ta-Tb)*vc;
    % if(vc < 0 || vc > vc_max || abs(l-pt) > cons || (T-Ta-Tb) < -cons)
    if(vc < -cons || vc > vc_max + cons || abs(l-pt) > max(pt,1.0)*1e-10)
        error('wrong in s_scurve_cpt_vavc CASE 1')
    end
    return;
end

% CASE 2: vc < min(vb,va_range)
% va = va_below
% vc < min(va_below, vb)
T_va_below_to_vb  = s_acc_time(va_below, vb, a, j);
l_below = min(va_below, vb) * (T - T_va_below_to_vb) + T_va_below_to_vb*(vb+va_below)/2;
if(pt < l_below)
    % 如果用 vc 做未知数，用 newton-ranphson 方法，在vc 很接近v2的时候数值性能很差
    % 例如：
    va = va_below;
    vc_below = s_cpt_vc_below_by_va_vb_T(va,vb,T,a,j);

    vc = newton_raphson_binary_search(@(vc)(...
        s_acc_time(va,vc,a,j)*(va+vc)/2  ...
        + s_acc_time(vb,vc,a,j)*(vb+vc)/2  ...
        + (T-s_acc_time(va,vc,a,j)-s_acc_time(vb,vc,a,j))*vc ...  % 本行不同
        - pt) ...
        ,max(vc_below, 0),min(va,vb)...
        ,10*eps);
    Ta = s_acc_time(va,vc,a,j);
    Tb = s_acc_time(vb,vc,a,j);
    mode = 0;
    l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0 + (T-Ta-Tb)*vc;
    
%     va = va_below;
%     v1 = min(va,vb);
%     v2 = max(va,vb);
%     
%     % 事实上 newton-ranphson 内的方程应该为：
%     % vc = s_acc_vend(v1,-a,-j,T1);
%     % T2 = s_acc_time(v2,vc,a,j);
%     %     
%     % T1*(v1+vc)/2 + T2*(v2+vc)/2 + (T-T2-T1)*vc - pt
%     
%     % 但matlab似乎只支持单行函数
%     T1_upper = min(newton_raphson_binary_search(@(T1)(...
%         s_acc_vend(v1,-a,-j,T1) - s_acc_vend(v2,-a,-j,T-T1)) ...
%         ,0,T...
%         ,eps), s_acc_time(v1,0,a,j));
%     
%     T1 = newton_raphson_binary_search(@(T1)(...
%         T1*(v1+s_acc_vend(v1,-a,-j,T1))/2  ...
%         + s_acc_time(v2,s_acc_vend(v1,-a,-j,T1),a,j)*(v2+s_acc_vend(v1,-a,-j,T1))/2  ...
%         + (T-T1-s_acc_time(v2,s_acc_vend(v1,-a,-j,T1),a,j))*s_acc_vend(v1,-a,-j,T1) ...  % 本行不同
%         - pt) ...
%         ,0,T1_upper...
%         ,10*eps);
% 
%     vc = s_acc_vend(v1,-a,-j,T1);
%     T2 = s_acc_time(v2,vc,a,j);
%     if(va < vb)
%         Ta = T1;
%         Tb = T2;
%     else
%         Tb = T1;
%         Ta = T2;
%     end
    

    % debug check %
    l = Ta* (va + vc) / 2.0 + Tb * (vb + vc) / 2.0 + (T-Ta-Tb)*vc;
    if(vc < 0 || vc > vc_max + cons || abs(l-pt) > 1e-10*max(pt,1))
        T1 = T1_upper;
        vc = s_acc_vend(v1,-a,-j,T1_upper);
        T2 = T-T1;
        l = T1* (v1 + vc) / 2.0 + T2 * (v2 + vc) / 2.0 + (T-T1-T2)*vc;
        l-pt
        T1 = T1_upper;
        vc = s_acc_vend(v1,-a,-j,T1_upper);
        T2 = s_acc_time(v2,vc,a,j);
        l = T1* (v1 + vc) / 2.0 + T2 * (v2 + vc) / 2.0 + (T-T1-T2)*vc;
        l-pt
        error('wrong in s_scurve_cpt_vavc CASE 2')
    end

    return;
end

% CASE 3     : min(vb,va_range) < vc < max(vb,va_range)
%      3.1   : vb < va_below
%      3.1.1 : vb < vc < va_below   
%      3.1.2 : va_below < vc < va_upper
%
%      3.2   : va_upper < vb
%      3.1.1 : va_upper < vc < vb   
%      3.1.2 : va_below < vc < va_upper
%
%      3.3   : va_below < vb < va_upper

% CASE 3.1
if(vb <= va_below)
    l_mid = T_va_below_to_vb*(va_below + vb)/2 + (T-T_va_below_to_vb)*va_below;
    if(pt < l_mid)
        % 3.1.1
        va = va_below;
        if(T-T_va_below_to_vb > 1e-10)
            vc = (pt - T_va_below_to_vb*(va + vb)/2)/(T-T_va_below_to_vb);
            Ta=abs((vc - vb)/(va - vb)) * (T-T_va_below_to_vb);
            Tb=abs((vc - va)/(va - vb)) * (T-T_va_below_to_vb);
        else
            vc = (vb+va_below)/2;
            Ta = (T-T_va_below_to_vb)/2;
            Tb = (T-T_va_below_to_vb)/2;
        end
        mode = 1;
        return;
    else
        % 3.1.2
        % va * (T-Tb) + (va+vb)/2*Tb == pt
        %
        va = newton_raphson_binary_search(@(va)(...
            s_acc_time(vb,va,a,j)*(va+vb)/2  ...
            + (T-s_acc_time(va,vb,a,j))*va ...  
            - pt) ...
            ,va_below,va_upper...
            ,10*eps);
        vc = va;
        Ta = T - s_acc_time(vb,va,a,j);
        Tb = 0;
        mode = 1;
        return;
    end
elseif(vb >= va_upper)
    l_mid = T_va_upper_to_vb*(va_upper + vb)/2 + (T-T_va_upper_to_vb)*va_upper;
    if(pt > l_mid)
        % 3.2.1
        va = va_upper;
        if(T-T_va_upper_to_vb > 1e-10)
            vc = (pt - T_va_upper_to_vb*(va_upper + vb)/2)/(T-T_va_upper_to_vb);
            Ta=abs((vc - vb)/(va - vb)) * (T-T_va_upper_to_vb);
            Tb=abs((vc - va)/(va - vb)) * (T-T_va_upper_to_vb);
        else
            vc = (vb+va_upper)/2;
            Ta = (T-T_va_upper_to_vb)/2;
            Tb = (T-T_va_upper_to_vb)/2;
        end
        mode = 1;
        return;
    else
        % 3.2.2
        % va * (T-Tb) + (va+vb)/2*Tb == pt
        %
        va = newton_raphson_binary_search(@(va)(...
            s_acc_time(vb,va,a,j)*(va+vb)/2  ...
            + (T-s_acc_time(va,vb,a,j))*va ... 
            - pt) ...
            ,va_below,va_upper...
            ,10*eps);
        vc = va;
        Ta = T - s_acc_time(vb,va,a,j);
        Tb = 0;
        mode = 1;
        return;
    end
else
    va = newton_raphson_binary_search(@(va)(...
        s_acc_time(vb,va,a,j)*(va+vb)/2  ...
        + (T-s_acc_time(va,vb,a,j))*va - pt), ...
        va_below,va_upper, ...
        10*eps);
    vc = va;
    Ta = T - s_acc_time(vb,va,a,j);
    Tb = 0;
    mode = 1;
    return;
end


end


function vc = s_scurve_cpt_vavc_case1(pt, va, vb, vc_max, a, j, T)
    % CASE 1.1 : va -> vc 无匀速段，vc -> vb 无匀速段
    %      1.2 : va -> vc 无匀速段，vc -> vb 有匀速段
    %      1.3 : va -> vc 有匀速段，vc -> vb 无匀速段
    %      1.4 : va -> vc 有匀速段，vc -> vb 有匀速段
    
    vc_upper = s_cpt_vc_upper_by_va_vb_T(va,vb,T,a,j);

    vc = newton_raphson_binary_search(@(vc)(...
        s_acc_time(va,vc,a,j)*(va+vc)/2  ...
        + s_acc_time(vb,vc,a,j)*(vb+vc)/2  ...
        + max(T-s_acc_time(va,vc,a,j)-s_acc_time(vb,vc,a,j), 0)*vc ...
        - pt) ...
        ,max(va,vb),min(vc_max, vc_upper)...
        ,10*eps);
    return;
end

function vc = s_scurve_cpt_vavc_case2(pt, va, vb, vc_max, a, j, T)
    % CASE 1.1 : va -> vc 无匀速段，vc -> vb 无匀速段
    %      1.2 : va -> vc 无匀速段，vc -> vb 有匀速段
    %      1.3 : va -> vc 有匀速段，vc -> vb 无匀速段
    %      1.4 : va -> vc 有匀速段，vc -> vb 有匀速段
    
    vc_below = s_cpt_vc_below_by_va_vb_T(va,vb,T,a,j);

    vc = newton_raphson_binary_search(@(vc)(...
        s_acc_time(va,vc,a,j)*(va+vc)/2  ...
        + s_acc_time(vb,vc,a,j)*(vb+vc)/2  ...
        + (T-s_acc_time(va,vc,a,j)-s_acc_time(vb,vc,a,j))*vc ...  % 本行不同
        - pt) ...
        ,max(vc_below, 0),min(va,vb)...
        ,10*eps);
    return;
end

