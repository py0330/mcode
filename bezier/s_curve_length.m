function [p,v_real] = s_curve_length(va, vb, v, a, j, T)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

% va : 起始速度
% vb : 终止速度
% v  : 最大匀速度（可能达不到），一定有 v >= vb
% a  : 最大加速度
% j  : 最大加加速度
% T  : 时间长度
%
% 本函数严格保证以下条件：
% 1. 起始速度为 va，结束速度为 vb
% 2. 起始加速度为0，结束加速度为0
% 3. 过程中最大加速度不超过 a
% 4. 过程种最大加加速度为 j
% 5. 所花时间为 T
% 6. v >= vb
%
% 本函数不保证如下条件：
% 1. 过程中一定可以达到最大速度v
%
% 函数返回：
% 实际的轨迹长度 : p
% 实际的最大速度 : v_real


% 计算满足va vb的曲线长度
% 需分2种情况考虑：
% ----------------------------------------------------
% CASE 1：v > max(va,vb)
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
% CASE 2：v在va,vb之间，例如va > v > vb
%
% va -------
%             \       
%                \          
%                   \        
%                      ------- vb
%  |<- Tb ->|<- Ta ->|<- Tc ->|
% 此时，轨迹仅有1段降速过程，即从va降速到vb，
% 但在va和vb处，会各有一段匀速过程
% 这两段匀速过程持续的时间 按va-v 和 vb-v的比例来决定

cons = eps * 10;

% 确定根据起始终止速度要求，可以满足时间T
if(s_acc_time(va,vb,a,j) > T + 1e-10)
    % ERROR! FAILED BECAUSE TIME INVALID
    % 需确保在T时间内可以从 va 加（减）速到vb
    p=-1;
    v_real = -1;
    return;
end


if(v >= max(va,vb))
    % CASE 1
    % 确定是否可以达到指定的v
    v1 = max(va,vb);
    v2 = min(va,vb);

    % 总体过程应分三种情况：
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
    % v  = v2+a^2/j
    % t1 = 2 * sqrt((v-v1) / j )
    % p1 = (v+v1)/2*t1
    % t2 = (v-v2) / a + a/j
    % p2 = (v+v2)/2*t2
    % T_ = t1 + t2
    %    = 2*((a^2/j - v1 + v2)/j)^(1/2) + (2*a)/j
    % 当 v1 > v2+a^2/j时，或 T > T_时
    % 若v1 < v2+a^2/j 且 T > T_，则处于case2
    %
    % CASE 3:
    % 其余情况
    if(s_acc_time(va,v,a,j) + s_acc_time(vb,v,a,j) > T)
        if(T > (4*a/j + v1/a - v2/a))
            % 此时有：
            v = min((T - (2*a)/j + v1/a + v2/a)/2*a, v);
        elseif(v1 > v2+a^2/j || T > (2*((a^2/j - v1 + v2)/j)^(1/2) + (2*a)/j))
            % 此时：
            % t1 = 2 * sqrt((v-v1) / j )
            % t2 = (v-v2) / a + a/j
            % t1+t2=T
            % 因此列方程
            % t1^2 - (T-t2)^2 = 0
            % 为关于v的一元二次方程：
            % 【1,- 2*v2 - 2*T*a - (2*a^2)/j,a^2*((T - a/j + v2/a)^2 + (4*v1)/j)】
            ratio_b = - 2*v2 - 2*T*a - (2*a^2)/j;
            ratio_c = a^2*((T - a/j + v2/a)^2 + (4*v1)/j);
            v_solution1 = (-ratio_b - sqrt(ratio_b*ratio_b - 4*ratio_c))/2;
            v = min(v_solution1,v);
        else
            % clear
            % syms va j pt T v vb f(Ta)
            % Ta = 2 * sqrt((v-va) / j )
            % Tb = 2 * sqrt((v - vb) / j )
            % eq = (Ta+Tb)^2
            % left = 8*(v/j - va/j)^(1/2)*(v/j - vb/j)^(1/2)
            % right = T^2 + (4*va)/j + (4*vb)/j - (8*v)/j
            % expand(left^2 - right^2)
            % collect(left^2 - right^2,v)
            % expand(((16*(T^2 + (4*va)/j + (4*vb)/j))/j - (64*va)/j^2 - (64*vb)/j^2))
            % expand((64*va*vb)/j^2 - (T^2 + (4*va)/j + (4*vb)/j)^2)
            v = - T^4*j^2 - 8*T^2*j*v1 - 8*T^2*j*v2 - 16*v1^2 + 32*v1*v2 - 16*v2^2;
            v = -v/(16*T^2*j);
        end
    end

    Ta = s_acc_time(va,v,a,j);
    Tb = s_acc_time(vb,v,a,j);

    p = (va + v)/2 * Ta + (vb + v)/2 * Tb + v * (T - Ta - Tb);
    v_real =v;
else
    % CASE 2
    Tc = s_acc_time(va,vb,a,j);
    total_even_time = T - Tc;
    Ta = (1-abs(va-v)/abs(va-vb)) * total_even_time;
    Tb = (1-abs(vb-v)/abs(va-vb)) * total_even_time;
    
    p = va * Ta + vb * Tb + (va + vb)/2 * Tc;
    v_real = v;
end



end