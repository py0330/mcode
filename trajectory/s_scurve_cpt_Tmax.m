function Tmax = s_scurve_cpt_Tmax(pa, va, pb, vb_max, vc_max, a, j)
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
if(va > vb_max && l_va_to_vb > pt + cons)
    Tmax=-1;
    Tmin=-1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 计算Tmax：%%%%%%%%
%
% 计算在最大的减速到0的情况下，所可能达到的最大长度
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



% 计算vb
% 因为可能有多解需要计算出距离 vb_max 最近的 vb 值
%
% 1）在 vb >= va - a^2/j 时，此时没有匀加速段
% 此时无匀加速段
% 前进时间为：t = 2 * sqrt( (va-vb) / j );
% 前进长度的平方为：
% p^2 = (t * (va+vb)/2)^2
%     =  va^3/j + (va^2*vb)/j - (va*vb^2)/j - vb^3/j
% vb 需求解一元三次方程 【1,va,-va^2,pb-pa-va^3】
%
% vb 取尽可能大的实数
% vb 的范围取自 【va/3，va】
%
% 1.1）va/3 <= va - a^2/j
% l 随 vb 在[va/3，va]单调递减，极大值位于 vb = va/3 处，极小值位于 vb = va处
% 考虑 vb 应当尽量贴近 vb_max，因此不考虑 vb < va/3 时的根
%
% 1.2）va/3 > va - a^2/j
% l 随 vb 在[va - a^2/j，va]单调递减，极大值位于 vb = va - a^2/j 处，极小值位于 vb = va处
%
% 2）在 vb <  va - a^2/j 时，此时拥有匀加速段
% 此时有匀速段
% 前进时间为：t = (va-vb)/a+a/j;
% 前进长度为：
% l = t*(va+vb)/2
%   = - vb^2/(2*a) + (a*vb)/(2*j) + (va*(a/j + va/a))/2
% vb 需求解一元二次方程 【
%       1,
%       (va - a*(a/j + va/a)),
%       2*(pb-pa)*a- a*va*(a/j + va/a)
% 】
% 对于根来说，应当取大值，这是因为Tmax应该尽可能的小
% 其极值应当位于 a^2/(2*j) 处
%
% 2.1）a^2/(2*j) <= va - a^2/j
% l 随 vb 在[a^2/(2*j)，va - a^2/j]单调递减，极大值位于 vb = a^2/(2*j) 处，极小值位于 vb = va - a^2/j 处
%            考虑 vb 应当尽量贴近 vb_max，因此不考虑 vb < a^2/(2*j) 时的根
%
% 2.2）a^2/(2*j) >  va - a^2/j
% l 随 vb 在[0, va - a^2/j]单调递增，极大值位于 vb = va - a^2/j 处，极小值位于 vb = 0处

% 【条件1】 加速度正好可以达到a时，所前进的长度
% 此时 vb = - a^2/j + va
% 前进时间 t = (va-vb)/a+a/j
% 前进长度 l = t*(va+vb)/2 = (2*a*va)/j - a*a*a/j/j

vb = max(va/3,va-a^2/j);
if(vb <= vb_max)
    lmax1       = s_acc_time(va,vb,a,j) * (va + vb)/2;
    if(lmax1>pt)
        vb = newton_raphson_binary_search(@(x)(sqrt((va-x)/j) * (va+x) - pt)...
                ,vb,va,10*eps);
        Tmax = s_acc_time(va,vb,a,j);
        return;
    end
end

vb = min([a^2/(2*j),va-a^2/j,vb_max]);
if(vb >= 0 && vb <= vb_max)
    lmax2       = s_acc_time(va,vb,a,j) * (va + vb)/2;
    if(lmax2>pt)
        B = -Z1;
        C = 2*pt*a- va*Z1 - va*va;
        vb = (-B + sqrt(B*B-4*C))/2;
        Tmax = s_acc_time(va,vb,a,j);
        return;
    end
end

Tmax = inf;
return;


if(va - 3/2*Z1 >  0)
    pacc = 1/2*(Z1/2 + va)^2/a;
else
    pacc = 4/3*va*sqrt(2/3*va/j);
end

if(pacc <= pt)
    Tmax = inf;
else
    % 计算vb
    % 【条件1】 加速度正好可以达到a时，所前进的长度
    % 此时 vb = - a^2/j + va
    % 前进时间 t = (va-vb)/a+a/j
    % 前进长度 l = t*(va+vb)/2 = (2*a*va)/j - a*a*a/j/j
    if(va<Z1 || (2*a*va)/j - a*a*a/j/j > pt)
        % 此时无匀速段
        % 前进时间为：t = 2 * sqrt( (va-vb) / j );
        % 前进长度的平方为：
        % p^2 = (t * (va+vb)/2)^2
        %     =  va^3/j + (va^2*vb)/j - (va*vb^2)/j - vb^3/j
        % vb 需求解一元三次方程 【1,va,-va^2,pb-pa-va^3】
        %
        % vb 取尽可能大的实数
        % vb 的范围取自 【va/3，va】,
        % 因为T的极值为 sqrt(va*8/3/j)，此时带入vb的公式，可得

        vb = newton_raphson_binary_search(@(x)(sqrt((va-x)/j) * (va+x) - pt)...
            ,va/3,va,10*eps);

        Tmax = s_acc_time(va,vb,a,j);
    else
        % 此时有匀速段
        % 前进时间为：t = (va-vb)/a+a/j;
        % 前进长度为：
        % l = t*(va+vb)/2
        %   = - vb^2/(2*a) + (a*vb)/(2*j) + (va*(a/j + va/a))/2
        % vb 需求解一元二次方程 【
        %       1,
        %       (va - a*(a/j + va/a)),
        %       2*(pb-pa)*a- a*va*(a/j + va/a)
        % 】
        % 对于根来说，应当取大值，这是因为Tmax应该尽可能的小
        B = -Z1;
        C = 2*pt*a- va*Z1 - va*va;
        vb = (-B + sqrt(B*B-4*C))/2;
        Tmax = s_acc_time(va,vb,a,j);
    end
end

end
