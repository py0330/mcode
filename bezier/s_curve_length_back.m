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