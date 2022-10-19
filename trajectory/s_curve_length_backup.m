function [p,v_real] = s_curve_length(va, vb, v, a, j, T)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

% 计算满足va vb的曲线长度
% 需分2种情况考虑：
% ----------------------------------------------------
% CASE 1：v > max(va,vb) or v < min(va,vb) 
%
%          ---                 v
%        /     \
% va ---        \
%                \
%                  --- vb      
%
% 此时，真实的v可能会无法取到，因而real_v可能会改变. 
% 若 v > max(va,vb), 则real_v 在会取[max(va,vb) v] 中的最大值
% 若 v < min(va,vb), 则real_v 在会取[min(va,vb) v] 中的最小值
% -----------------------------------------------------
% CASE 2：v在va,vb之间，例如va > v > vb
%
% va ---
%        \       
%          ---                 v  
%              \
%                \
%                  --- vb
%
% 此时，轨迹仅有1段降速过程，即从va降速到vb，但在va和vb处，会有一段匀速过程4
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


if(v >= max(va,vb) || v < min(va,vb))
    % CASE 1
    % 确定是否可以达到指定的v

    if(v >= max(va,vb))
        v1 = min(va,vb);
        v2 = max(va,vb);
        
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
        if(s_acc_time(va,v,a,j) + s_acc_time(vb,v,a,j) > T)
        if(T > (4*a/j + v1/a - v2/a))
            % 此时有：
            v = min((T - (2*a)/j + v1/a + v2/a)/2*a, v);
%         elseif(v1 < v2+a^2/j && T > (2*((a^2/j - v1 + v2)/j)^(1/2) + (2*a)/j))
            % 此时：
            % t1 = 2 * sqrt((v-v1) / j )
            % t2 = (v-v2) / a + a/j
            % t1+t2=T
            % 因此列方程
            % t1^2 - (T-t2)^2 = 0
            % 为关于v的一元二次方程：
            % 【1,- 2*v2 - 2*T*a - (2*a^2)/j,a^2*((T - a/j + v2/a)^2 + (4*v1)/j)】
%             ratio_b = - 2*v2 - 2*T*a - (2*a^2)/j;
%             ratio_c = a^2*((T - a/j + v2/a)^2 + (4*v1)/j);
%             v = (-ratio_b + sqrt(ratio_b*ratio_b - 4*ratio_c))/2;
        else
            v
            if(s_acc_time(va,v,a,j) + s_acc_time(vb,v,a,j) > T)
                if(v > (va+vb)/2)
                    v_success = max(va,vb);
                    v_failed  = v;
                else
                    v_success = min(va,vb);
                    v_failed  = v;
                end

                while(abs(v_success - v_failed)>cons)
                    v_next = (v_success + v_failed)/2;
                    if(s_acc_time(va,v_next,a,j) + s_acc_time(vb,v_next,a,j) > T)
                        v_failed  = v_next;
                    else
                        v_success = v_next;
                    end
                end
                v = v_success;
            end

        end
        end
    else
        if(s_acc_time(va,v,a,j) + s_acc_time(vb,v,a,j) > T)
            if(v > (va+vb)/2)
                v_success = max(va,vb);
                v_failed  = v;
            else
                v_success = min(va,vb);
                v_failed  = v;
            end

            while(abs(v_success - v_failed)>cons)
                v_next = (v_success + v_failed)/2;
                if(s_acc_time(va,v_next,a,j) + s_acc_time(vb,v_next,a,j) > T)
                    v_failed  = v_next;
                else
                    v_success = v_next;
                end
            end
            v = v_success;
        end
    end
    Ta = s_acc_time(va,v,a,j);
    Tb = s_acc_time(vb,v,a,j);

    p = (va + v)/2 * Ta + (vb + v)/2 * Tb + v * (T - Ta - Tb);
    v_real =v;
else
    % CASE 2
    total_even_time = T - s_acc_time(va,vb,a,j);
    Ta = (1-abs(va-v)/abs(va-vb)) * total_even_time;
    Tb = (1-abs(vb-v)/abs(va-vb)) * total_even_time;
    
    p = va * Ta + vb * Tb + (va + vb)/2 * (T - Ta - Tb);
    v_real = v;
end



end