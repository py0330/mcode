function [vc_upper,Ta,Tb] = s_cpt_vc_upper_by_va_vb_T(va,vb,T,a,j)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
    
    cons = 1e-10;

    v1 = min(va,vb);
    v2 = max(va,vb);
    ve1 = v1 + a*a/j;
    ve2 = v2 + a*a/j;

    T_v1_to_v2 = s_acc_time(v1,v2,a,j);
    T_v1_to_ve1 = s_acc_time(v1,ve1,a,j);
    T_v2_to_ve1 = s_acc_time(v2,ve1,a,j);
    T_v1_to_ve2 = s_acc_time(v1,ve2,a,j);
    T_v2_to_ve2 = s_acc_time(v2,ve2,a,j);

    % v1 无法加速到 v2
    if(T_v1_to_v2 > T + cons)
        error('failed in s_cpt_vc_upper_by_va_vb_T')
    end

    if(ve1 >= v2 && T_v1_to_ve1 + T_v2_to_ve1 >= T)
        % 第一段无匀加速，第二段无匀加速
        % clear
        % syms T1 T2 T v1 v2 vc a j pt
        % T1 = 2 * sqrt((vc-v1)/j)
        % T2 = 2 * sqrt((vc-v2)/j)
        % vc_ans = solve(T1 + T2 == T, vc)
        vc = (T^4*j^2 + 8*T^2*j*v1 + 8*T^2*j*v2 + 16*v1^2 - 32*v1*v2 + 16*v2^2)/(16*T^2*j);
        T1 = 2 * sqrt((vc-v1)/j);
        T2 = 2 * sqrt((vc-v2)/j);
    elseif(T_v1_to_ve2 + T_v2_to_ve2 >= T)
        % 第一段有匀加速，第二段无匀加速
        % clear
        % syms T1 T2 T v1 v2 vc a j pt
        % T1 = (vc-v1)/a + a/j
        % T2 = 2 * sqrt((vc-v2)/j)
        % solve(T1 + T2 == T, vc)
        vc = (j*v1 + a^2 - 2*a*j*((v1 - v2 + T*a)/j)^(1/2) + T*a*j)/j;
        T1 = (vc-v1)/a + a/j;
        T2 = 2 * sqrt((vc-v2)/j);
    else
        % 第一段有匀加速，第二段有匀加速
        % clear
        % syms T1 T2 T v1 v2 vc a j pt
        % T1 = (vc-v1)/a + a/j
        % T2 = (vc-v2)/a + a/j
        % solve(T1 + T2 == T, vc)
        vc = (- 2*a^2 + T*j*a + j*v1 + j*v2)/(2*j);
        T1 = (vc-v1)/a + a/j;
        T2 = (vc-v2)/a + a/j;
    end
    
    
    % 更新 vc_upper, Ta, Tb
    vc_upper = vc;
    if(vb > va)
        Ta = T1;
        Tb = T2;
    else
        Ta = T2;
        Tb = T1;
    end
end
