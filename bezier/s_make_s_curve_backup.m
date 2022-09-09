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

% 二分法找vb，尽可能让结束时的速度最小
if(s_curve_length(va, vb_min, max_v, a, j, T) > pt)
    vb = vb_min;
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

% s_curve_vrange(va, vb, a, j, T)

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



