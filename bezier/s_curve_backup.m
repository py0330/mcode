function [p, v1] = s_curve(t, p0, v0, p1, a, j, T)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
% init acc: zero
%
% t      : current time
% p0     : init pos
% v0     : init vel
% p1     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period

% si 如下：si = sign(p1 - T*v0)     |   若（p1 - T*v0 > 0）则 v1 > v0
pt = p1 - p0;
si = sign(pt - T*v0);

% 应该有：p1 = p0 + (v0+v1)/2 * T1 + v1 * T2           【等式1】
% T1 为加速段时间，将速度从v0 加速到 v1
% T2 为匀速段时间，以v1速度匀速运动到结束
% 因此：T1 + T2 = T

% 下推 v1 与 T1 的关系

% CASE 1 下，T1的计算方案如下：
% CASE 1 |v1 - v0| > a^2 / j, 此时加速段可以达到最大加速度：
%        T1 = |v0 - v1|/a + a/j            |    |v1 - v0| >  a^2 / j
%           = si/a*(v1 - v0) + a/j
%        v1 = si*a*T1 + (v0 - si*a^2/j)        【等式3】
%
%        将上式代入【等式1】中，可得：
%   a/2 * T1^2  -(a^2/(2*j) + 3*a) * T1 + si*(p1 - p0 - 3*v0) + (3*a^2)/j == 0
%
%        求取上述一元二次方程，可得T1


% T1_candidate = square_equation_solve(1, -a/j - 2*T, 2*si*(pt - T*v0)/a + 2*T*a/j);
% for i =1:2
%     if(T1_candidate(i) > T || T1_candidate(i) < 0)
%         continue
%     else
%         T1 = T1_candidate(i);
%     end
% end


square_a = 1;
square_b = -a/j - 2*T;
square_c = 2*si*(pt - T*v0)/a + 2*T*a/j;

T1_candidate = (-square_b - sqrt(square_b*square_b-4*square_a*square_c))/(square_a*2);

if(isreal(T1_candidate) && T1_candidate >= 2*a/j && T1_candidate <= T)
    T1 = T1_candidate;

    v1 = si*a*T1 + (v0 - si*a^2/j);
    if(t < a/j)
        p = p0 + v0 * t + si/6*j*t^3;
    elseif(t < (T1 - a/j))
        p = p0 + v0 * a / j + si/6 * a^3 / j^2 + (v0 + si/2*a^2/j) * (t - a/j)+ si*a/2*(t-a/j)^2; 
    elseif(t < T1)
        p = p0 + (v0+v1)/2*T1-(v1*(T1-t) - si/6*j*(T1-t)^3);
    else
        p = p0 + (v0+v1)/2*T1 + v1 * (t - T1);
    end
    
    return;
end

% CASE 2 若 |v1 - v0| <= a^2 / j, 此时加速段（或减速段）无法达到最大加速度：
%        T1 = 2 * sqrt( |v0 - v1| / j )    |    |v1 - v0| <= a^2 / j
%        因此：
%        v1 = si/4 * j * T1^2 + v0                     【等式2】
%        
%        将 v1 带入到【等式1】中：
%        T1^3 - 2*T* T1^2 + si*8/j*(pt - T*v0) == 0
%        解该一元三次方程，可得T1
%        合理的T1应该是： 0 <= T1 <= T 且 j*T1/2 < a
%        此后可以计算出 v1
% CASE 2 下，T1的计算方案如下：
T1_candidate = cubic_equation_solve(1, -2*T, 0, si*8/j*(pt - T*v0));
T1_candidate = sort(T1_candidate);

% 因为在 T1 = 0 的时候，原等式> 0, 因此一元三次方程必有一个负根。
% 因此如果方程有正根，那么一定有2个正根，T1一定是小的那个
T1 = T1_candidate(2);

v1 = si/4 * j * T1^2 + v0;

if(t < T1 / 2)
    p = p0 + v0 * t + si/6*j*t^3;
elseif(t < T1)
    p = p0 + (v0+v1)/2*T1 - (v1 * (T1-t) - si/6*j*(T1-t)^3);
else
    p = p0 + (v0+v1)/2*T1 + v1 * (t - T1);
end

end



