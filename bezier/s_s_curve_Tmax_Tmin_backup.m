function [Tmax, Tmin] = s_s_curve_Tmax_Tmin(pa, va, pb, max_vb, max_v, a, j)
% 计算当前点位所需的最大最小时间
%
% pa     : init pos
% va     : init vel
% pb     : end  pos
% max_vb : max end vel
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period
%
% Tmax：开始时尽可能快的减速，若减速到0，则为inf，否则以到达pb的时间为准
% Tmin：开始时尽可能快的加速，直到速度最大，之后保持最大速度到终点

cons = eps * 10;

% 计算Tmax：
%
% 计算在最大的减速情况下，所可能达到的最大长度
% 【注】：最大长度并非发生在最大的T处
% 即：p并非是T的单调函数
%
% 下求最大的长度：
% 终止速度vb为：
% vb = va + a^2/j - a * T                when va - vb >  a^2 / j
%      va - j * T^2 / 4                  when va - vb <= a^2 / j
%      
% 此时行进的长度为：
% p = (va + vb)/2 * T
%   = -a/2 * T^2 + (va + a^2/j/2) * T    when va - vb >  a^2 / j
%     -j/8 * T^3 + va * T                when va - vb <= a^2 / j
%
% 对其求导，并求极值
% dp = -a*T + (va + a^2/j/2)             when va - vb >  a^2 / j
%      -3*j/8 * T^2 + va                 when va - vb <= a^2 / j
%
% 即：
% T = (a^2/(2*j) + va)/a                 when va - vb >  a^2 / j
%     sqrt(va*8/3/j)                     when va - vb <= a^2 / j
%
% 判断极值的触发条件：
% 将上述T，带入到右侧判别式，即 va - vb - a^2 / j 中，可以发现：
% va - vb - a^2 / j =     va - 3/2*a^2/j
% or
% va - vb - a^2 / j = 2/3*va -     a^2/j
%
% 两者仅相差一个系数，因此可能产生最长路径的T为：
% T = (a^2/(2*j) + va)/a                 when va - 3/2*a^2/j >  0
%     sqrt(va*8/3/j)                     when va - 3/2*a^2/j <= 0
%
% 带入可得此时的p：
% p = 1/2*(a^2/(2*j) + va))^2/a          when va - 3/2*a^2/j >  0
%   = 4/3*va*sqrt(2/3*va/j)              when va - 3/2*a^2/j <= 0

if(va - 3/2*a^2/j >  0)
    Tacc = (a^2/(2*j) + va)/a;
    pacc = 1/2*(a^2/(2*j) + va)^2/a;
else
    Tacc = sqrt(va*8/3/j);
    pacc = 4/3*va*sqrt(2/3*va/j);
end

% Tacc = s_acc_time(va,0,a,j);
% pacc = Tacc * (va + 0)/2;



if(pacc <= pb - pa + 1e-10)
    Tmax = inf;
else
    % 计算vb
    % 以下条件判断正好没有匀速段时，所前进的长度
    % 此时 vb = - a^2/j + va
    % 前进时间 t = (va-vb)/a+a/j
    % 前进长度 p = t*(va+vb)/2 = (2*a*va)/j - a*a*a/j/j
    if((2*a*va)/j - a*a*a/j/j > pb - pa)
        % 此时无匀速段
        % 前进时间为：t = 2 * sqrt( (va-vb) / j );
        % 前进长度的平方为：
        % p^2 = (t * (va+vb)/2)^2
        %     =  va^3/j + (va^2*vb)/j - (va*vb^2)/j - vb^3/j
        % vb 需求解一元三次方程 【1,va,-va^2,pb-pa-va^3】
        r = cubic_equation_solve(1,va,-va^2,(pb-pa)*(pb-pa)*j-va^3);
        vb = r(3);
        Tmax = s_acc_time(va,vb,a,j);
    else
        % 此时有匀速段
        % 前进时间为：t = (va-vb)/a+a/j;
        % 前进长度为：
        % p = t*(va+vb)/2
        %   = - vb^2/(2*a) + (a*vb)/(2*j) + (va*(a/j + va/a))/2
        % vb 需求解一元二次方程 【
        %       1,
        %       (va - a*(a/j + va/a)),
        %       2*(pb-pa)*a- a*va*(a/j + va/a)
        % 】
        % 对于根来说，应当取大值，这是因为Tmax应该尽可能的小
        ratio_b = (va - a*(a/j + va/a));
        ratio_c = 2*(pb-pa)*a- a*va*(a/j + va/a);
        vb = (-ratio_b + sqrt(ratio_b*ratio_b-4*ratio_c))/2;
        Tmax = s_acc_time(va,vb,a,j);
    end
end

% 计算Tmin：
Tacc = s_acc_time(va,max_vb,a,j);
pacc = Tacc * (va + max_vb)/2;

if(pacc > pb - pa) % vb无法达到max_vb
    if(va > max_vb)% 起始速度高于max_vb，末端速度无法满足要求
        Tmax=-1;
        Tmin=-1;
        return;
        error('failed in s_s_curve_Tmax_Tmin: can not stop to max_vb');
    else % 起始速度低于最大许可的终止速度，计算最大可能的vb
        % binary search
        vb_upper = max_vb;
        vb_below = va;
        
        while(abs(vb_upper - vb_below) > cons)
            vb_next = (vb_upper + vb_below)/2;
    
            if(s_acc_time(va,vb_next,a,j)*(vb_next + va)/2 < pb-pa)
                vb_below = vb_next;
            else
                vb_upper = vb_next;
            end
        end
        vb = vb_upper;% v取大值，这样可以让时间Tmin更大
        Tmin = s_acc_time(va,vb,a,j);
    end
else
    Ta = s_acc_time(va   ,max_v ,a,j);
    Tb = s_acc_time(max_v,max_vb,a,j);
    if(Ta*(va + max_v)/2 + Tb*(max_v+max_vb)/2 > pb-pa) % 无法达到max_v
        % binary search v
        v_upper = max_v;
        v_below = max(va,max_vb);
        while(abs(v_upper - v_below) > cons)
            v_next = (v_upper + v_below)/2;
    
            if(    s_acc_time(va,v_next,a,j)*(v_next + va)/2 + ...
                   s_acc_time(v_next,max_vb,a,j)*(v_next + max_vb)/2 ...
                   < pb-pa)
                v_below = v_next;
            else
                v_upper = v_next;
            end
        end
        v = v_upper; % v取大值，这样可以增加Tmin：因为va和vb都小于v

        Tmin = s_acc_time(va,v,a,j) + s_acc_time(max_vb,v,a,j);
    else%可以达到max_v
        Tmin = Ta + Tb + ...
            (pb- pa -Ta*(va + max_v)/2 - Tb*(max_v+max_vb)/2)/max_v;
    end


end

end

function t = s_acc_time(va,vb,a,j)
    % CASE 1 |v - va| >  a^2 / j, 此时加速段可以达到最大加速度：
    %        T0 = |v - va|/a + a/j            |    |v - va| >  a^2 / j
    % CASE 2 |v - va| <= a^2 / j, 此时加速段（或减速段）无法达到最大加速度：
    %        T1 = 2 * sqrt( |v - va| / j )    |    |v - va| <= a^2 / j

    if(abs(vb-va)>a^2 / j)
        t = abs(vb-va) / a + a/j;
    else
        t = 2 * sqrt( abs(vb-va) / j );
    end
end

