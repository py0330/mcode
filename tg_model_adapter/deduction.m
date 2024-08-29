%%
% 构造3次函数，使得起始、终止点处的位置、速度与需求一致
% p(t)   = A*t^3 + B*t^2 + C*t + D
% dp(t)  = 3*A*t^2 + 2*B*t + C
% d2p(t) = 6*A*t + 2*B
% 
% [                1 ] * [A] = [pb]
% |             1    |   |B|   |vb|
% | T^3    T^2  T  1 |   |C|   |pe|
% [ 3*T^2  2*T  1    ]   [D]   [ve]

syms pb vb pe ve T vmax amax

M = [0,0,0,1
     0,0,1,0
     T^3, T^2, T,1
     3*T^2, 2*T,1,0];

b = [pb;vb;pe;ve];

x = M\b;

A = x(1);
B = x(2);
C = x(3);
D = x(4);

k2 = 3*A;
k1 = 2*B;
k0 = C;

l1 = 6*A;
l0 = 2*B;
%% 以下确定 dp 的极值点 是否在 (0,T)区间中，最终版
% m 为速度 dp 的极值点
m  = -k1/2/k2;

% 下求若存在速度极值时，极值不超过最大速度的条件
% 极值处速度vm为：
vm = 3*A*m^2 + 2*B*m + C;

% 应有 -vmax < vm < vmax
% 即：
% vm - vmax < 0
% vm + vmax > 0
eq1 = vm - vmax
eq2 = vm + vmax

% 设若其分子的四个根分别为：
% T1 T2 T3 T4
% 
% 可从大到小验证以上四个根的m 是否在 (0,T)中

% 还需有 -amax < a0 < amax
%        -amax < aT < amax
% 即：
%       -amax < 2*B < amax
%       -amax < 6*A*T + 2*B < amax
% 可得4个方程：
% 
% amax*T^2 - (2*vb + 4*ve)*T - (6*pb - 6*pe) > 0
% amax*T^2 + (2*vb + 4*ve)*T + (6*pb - 6*pe) > 0
% amax*T^2 + (4*vb + 2*ve)*T + (6*pb - 6*pe) > 0
% amax*T^2 - (4*vb + 2*ve)*T - (6*pb - 6*pe) > 0


% 
% 下证不可能出现 D1 < 0 && D2 < 0
% -D1 - D2 =  (A1-A2)*(vb+ve)
%         =  -6*vmax*(vb+ve)^2 
%         < 0
%
% 易证，在E1 > 0 && E2 > 0 时，必有 D1 > 0 && D2 > 0。
% 
% COND 1：D1 > 0 && D2 > 0
%    此时必定可以实现的 T 为 (max(T1,T4,T5,T6,T7), inf)
%
% COND 2：D1 < 0 && D2 > 0
%    则 eq1 的最大的两个根为 Tk1 < Tk2
%       eq2 的最大的两个根为 Tk3 < Tk4
%    若 Tk2 > max(T2,T3) && Tk2 > Tk4
%       则可取范围时 (Tk4,inf)
%    否则
%         可取范围是 (max(T2,T3),inf)
% COND 3：D1 > 0 && D2 < 0
%    则 eq1 的最大的两个根为 Tk1 < Tk2
%       eq2 的最大的两个根为 Tk3 < Tk4
%    若 Tk4 > max(T2,T3) && Tk4 > Tk2
%       则可取范围时 (Tk2,inf)
%    否则
%       可取范围是 (max(T2,T3),inf)




% 
% 上述不等式两个根为：
% T4,T5
% T6,T7
% 此时 T 可取范围是以下两者求交：
%   若 A1*(vb+ve) < 0，则为(T4,T5) 否则 (0,T4) 并 (T5,inf)
%   若 A2 > 0，则为(T6,T7) 否则 (0,T6) 并 (T7,inf)
%
% COND 2 T1 < 0
%
% 上述两个不等式分别可化为：
% A1 * T^2 + B1 * T + C1 > 0
% A2 * T^2 + B2 * T + C2 < 0
%
% 其中：
% A1 = -(- vb^2 - vb*ve - 3*vmax*vb - ve^2 - 3*vmax*ve)
% B1 = -(6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve - 6*pb*vmax + 6*pe*vmax)
% C1 = -(-9*pb^2 + 18*pb*pe - 9*pe^2)
%
% A2 = -(- vb^2 - vb*ve + 3*vmax*vb - ve^2 + 3*vmax*ve)
% B2 = -(6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve + 6*pb*vmax - 6*pe*vmax)
% C2 = -(-9*pb^2 + 18*pb*pe - 9*pe^2)
%
% 此时 T 可取范围是以下两者求交：
%   若 A1 < 0，则为(T4,T5) 否则 (0,T4) 并 (T5,inf)
%   若 A2 > 0，则为(T6,T7) 否则 (0,T6) 并 (T7,inf)
%
% 求交结果为：
%   若 A1 < 0 && A2 > 0: (max(T4,T6), min(T5,T7))
%   若 A1 < 0 && A2 < 0: (T4, min(T5,T6)) 并 (T7,T5)
%   若 A1 > 0 && A2 > 0: (T4, min(T5,T6)) 并 (T7,T5)
%   若 A1 > 0 && A2 < 0: (0,min(T4,T6))   并 (max(T5,T7),inf)
%
%
%【综上2】，肯定不存在极值点的条件为：
% T_max_candidate1 = max(0,T2,T3)
%
% T_max_candidate1 = max(0,T2,T3)
%
% T_max_candidate = (T1 > 0 && A1 < 0 && T5 > max(T2,T3)) ? : ;



%% 以下确定 dp 的极值点 是否在 (0,T)区间中
% m 为速度 dp 的极值点
m  = -k1/2/k2;

% 先判断 m 与 0 的关系
% m: (T*(3*pb - 3*pe + 2*T*vb + T*ve))/(3*(2*pb - 2*pe + T*vb + T*ve))
% 判断 m 与 0 的关系，此时只需要确定 m 的正负，故把上式的除法修改为乘法
eq1 = ((3*pb - 3*pe + 2*T*vb + T*ve))*(3*(2*pb - 2*pe + T*vb + T*ve));
% eq1 的符号与 m 相同
% eq1 = E1 * (T - T1) * (T - T2);
% 其中：
%     E1 = (vb + ve) * (2*vb + ve)
%     T1 = (2*pe-2*pb)/(vb+ve)
%     T2 = (3*pe-3*pb)/(2*vb+ve)
%
% 易知，E1 与 T1 * T2 的符号相同 
%
% 若 E1 < 0, 则 T > max(T1,T2) 或 T < min(T1,T2) 时，(0,T)区间内不存在速度极值
% 若 E1 > 0, 则 min(T1,T2) < T < max(T1,T2) 时，(0,T)区间内不存在速度极值




% 判断 m 与 T 的关系：此时类似上式进行判断
% 构造 T - m = T - (T*(3*pb - 3*pe + 2*T*vb + T*ve))/(3*(2*pb - 2*pe + T*vb + T*ve))
%            = T * (1 - ((3*pb - 3*pe + 2*T*vb + T*ve))/(3*(2*pb - 2*pe + T*vb + T*ve))
%            = T * (3*pb - 3*pe + T*vb + 2*T*ve)  / (3*(2*pb - 2*pe + T*vb + T*ve))
% 判断 T - m 的关系，此时只需要确定 m 的正负，故把上式的除法修改为乘法
eq2 = (3*pb - 3*pe + T*vb + 2*T*ve)*(3*(2*pb - 2*pe + T*vb + T*ve));
% eq2 的符号与 T - m 相同
% eq2 = E2 * (T - T1) * (T - T3);
% 其中：
%     E2 = (vb + ve) * (vb + 2*ve)
%     T1 = (2*pe-2*pb)/(vb+ve)
%     T3 = (3*pe-3*pb)/(vb+2*ve)
% 易知，E2 与 T1 * T3 的符号相同 
% 若 E2 < 0, 则 T > max(T1,T3) 或 T < min(T1,T3) 时，(0,T)区间内不存在速度极值
% 若 E2 > 0, 则 min(T1,T3) < T < max(T1,T3) 时，(0,T)区间内不存在速度极值

% 统一讨论，(0,T)区间内不存在极值的条件，
%
% Cond1 E1 > 0 && E2 > 0
% 此时 T1 、T2、 T3 同号，T2、T3 必然分布在 T1 的两侧
%      因此，无速度极值的区间为:
%      T1 > 0 时，(min(T1,T2,T3), max(T1,T2,T3))
%      T1 < 0 时，空集
%
% Cond2 E1 < 0 && E2 > 0
% 此时 T1 、T3 同号，T1、T2异号，T1、T2 必然分布在 T3 的两侧
%      因此，无速度极值的区间为:
%      T1 > 0 时，(T3,inf) 
%      T1 < 0 时，(T2,inf)
%
% Cond3 E1 > 0 && E2 < 0
% 此时 T1 、T2 同号，T1、T3异号，T1、T3 必然分布在 T2 的两侧
%      因此，无速度极值的区间为:
%      T1 > 0 时，(T2,inf) 
%      T1 < 0 时，(T3,inf)
%【综上1】，肯定不存在极值点的条件为：
% T 位于 (max(T2,T3), inf) 中

% 下求若存在速度极值时，极值不超过最大速度的条件
% 极值处速度vm为：
% vm = 3*A*m^2 + 2*B*m + C
vm = 3*A*m^2 + 2*B*m + C

% 应有 -vmax < vm < vmax
% 即：
% vm - vmax < 0
% vm + vmax > 0
eq1 = vm - vmax
eq2 = vm + vmax

% 上式分母处应为 (3*T*(2*pb - 2*pe + T*vb + T*ve))
% 即 (3*vb + 3*ve)*T^2 + (6*pb - 6*pe)*T
% 即 (3*vb + 3*ve)*T*(T-T1)
% 显然，其一个根为 0 ，另一个根应为2*(pe-pb)/(vb+ve)，即 T1
%
% 上述两个不等式分别可化为：
% (A1 * T^2 + B1 * T + C1)/((3*vb + 3*ve)*T*(T-T1)) < 0
% (A2 * T^2 + B2 * T + C2)/((3*vb + 3*ve)*T*(T-T1)) > 0
%
% 其中：
% A1 = (- vb^2 - vb*ve - 3*vmax*vb - ve^2 - 3*vmax*ve)
% B1 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve - 6*pb*vmax + 6*pe*vmax)
% C1 = - 9*pb^2 + 18*pb*pe - 9*pe^2
%
% A2 = (- vb^2 - vb*ve + 3*vmax*vb - ve^2 + 3*vmax*ve)
% B2 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve + 6*pb*vmax - 6*pe*vmax)
% C2 = - 9*pb^2 + 18*pb*pe - 9*pe^2
%
%
% 设若其分子的两个根分别为：
% T4, T5
% T6, T7
%
% 则上式可化为：
%
% D1 * (T-T4)*(T-T5)/(T*(T-T1)) > 0
% D2 * (T-T6)*(T-T7)/(T*(T-T1)) > 0
%
% 其中:
% D1 = -A1*(vb+ve)
% D2 = A2*(vb+ve)
%
% 进一步等同于
% D1 * T*(T-T1)*(T-T4)*(T-T5) > 0            eqn 1
% D2 * T*(T-T1)*(T-T6)*(T-T7) > 0            eqn 2
% 
% 下证不可能出现 D1 < 0 && D2 < 0
% -D1 - D2 =  (A1-A2)*(vb+ve)
%         =  -6*vmax*(vb+ve)^2 
%         < 0
%
% 易证，在E1 > 0 && E2 > 0 时，必有 D1 > 0 && D2 > 0。
% 
% COND 1：D1 > 0 && D2 > 0
%    此时必定可以实现的 T 为 (max(T1,T4,T5,T6,T7), inf)
%
% COND 2：D1 < 0 && D2 > 0
%    则 eq1 的最大的两个根为 Tk1 < Tk2
%       eq2 的最大的两个根为 Tk3 < Tk4
%    若 Tk2 > max(T2,T3) && Tk2 > Tk4
%       则可取范围时 (Tk4,inf)
%    否则
%         可取范围是 (max(T2,T3),inf)
% COND 3：D1 > 0 && D2 < 0
%    则 eq1 的最大的两个根为 Tk1 < Tk2
%       eq2 的最大的两个根为 Tk3 < Tk4
%    若 Tk4 > max(T2,T3) && Tk4 > Tk2
%       则可取范围时 (Tk2,inf)
%    否则
%       可取范围是 (max(T2,T3),inf)




% 
% 上述不等式两个根为：
% T4,T5
% T6,T7
% 此时 T 可取范围是以下两者求交：
%   若 A1*(vb+ve) < 0，则为(T4,T5) 否则 (0,T4) 并 (T5,inf)
%   若 A2 > 0，则为(T6,T7) 否则 (0,T6) 并 (T7,inf)
%
% COND 2 T1 < 0
%
% 上述两个不等式分别可化为：
% A1 * T^2 + B1 * T + C1 > 0
% A2 * T^2 + B2 * T + C2 < 0
%
% 其中：
% A1 = -(- vb^2 - vb*ve - 3*vmax*vb - ve^2 - 3*vmax*ve)
% B1 = -(6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve - 6*pb*vmax + 6*pe*vmax)
% C1 = -(-9*pb^2 + 18*pb*pe - 9*pe^2)
%
% A2 = -(- vb^2 - vb*ve + 3*vmax*vb - ve^2 + 3*vmax*ve)
% B2 = -(6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve + 6*pb*vmax - 6*pe*vmax)
% C2 = -(-9*pb^2 + 18*pb*pe - 9*pe^2)
%
% 此时 T 可取范围是以下两者求交：
%   若 A1 < 0，则为(T4,T5) 否则 (0,T4) 并 (T5,inf)
%   若 A2 > 0，则为(T6,T7) 否则 (0,T6) 并 (T7,inf)
%
% 求交结果为：
%   若 A1 < 0 && A2 > 0: (max(T4,T6), min(T5,T7))
%   若 A1 < 0 && A2 < 0: (T4, min(T5,T6)) 并 (T7,T5)
%   若 A1 > 0 && A2 > 0: (T4, min(T5,T6)) 并 (T7,T5)
%   若 A1 > 0 && A2 < 0: (0,min(T4,T6))   并 (max(T5,T7),inf)
%
%
%【综上2】，肯定不存在极值点的条件为：
% T_max_candidate1 = max(0,T2,T3)
%
% T_max_candidate1 = max(0,T2,T3)
%
% T_max_candidate = (T1 > 0 && A1 < 0 && T5 > max(T2,T3)) ? : ;

%%
% syms pb vb pe ve T vmax
vmax = 100;

pb = (rand-0.5) /(rand + 0.01);
pe = (rand-0.5) /(rand + 0.01);
vb = (rand-0.5) /(rand + 0.01);
ve = (rand-0.5) /(rand + 0.01);

E1 = (vb + ve) * (2*vb + ve);
E2 = (vb + ve) * (vb + 2*ve);

T1 = (2*pe-2*pb)/(vb+ve);
T2 = (3*pe-3*pb)/(2*vb+ve);
T3 = (3*pe-3*pb)/(vb+2*ve);

A1 = (- vb^2 - vb*ve - 3*vmax*vb - ve^2 - 3*vmax*ve);
B1 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve - 6*pb*vmax + 6*pe*vmax);
C1 = - 9*pb^2 + 18*pb*pe - 9*pe^2;

A2 = (- vb^2 - vb*ve + 3*vmax*vb - ve^2 + 3*vmax*ve);
B2 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve + 6*pb*vmax - 6*pe*vmax);
C2 = - 9*pb^2 + 18*pb*pe - 9*pe^2;

if(B1*B1-4*A1*C1 > 0)
    T_candidate1 = (-B1-sqrt(B1*B1-4*A1*C1))/(2*A1);
    T_candidate2 = (-B1-sqrt(B1*B1+4*A1*C1))/(2*A1);
    T4 = min([T_candidate2,T_candidate1]);
    T5 = max([T_candidate2,T_candidate1]);
else
    T4 = 0;
    T5 = inf;
end

if(B2*B2-4*A2*C2 > 0)
    T_candidate1 = (-B2-sqrt(B2*B2-4*A2*C2))/(2*A2);
    T_candidate2 = (-B2-sqrt(B2*B2+4*A2*C2))/(2*A2);
    T6 = min([T_candidate2,T_candidate1]);
    T7 = max([T_candidate2,T_candidate1]);
else
    T6 = 0;
    T7 = inf;
end





%%
eq = 3*A*m^2 + 2*B*m + C - vmax

eq1 = eq * (3*T*(2*pb - 2*pe + T*vb + T*ve));
expand(eq1)
collect(eq1,T)

collect(3*T*(2*pb - 2*pe + T*vb + T*ve),T)


%%
for i=1:1e4
    vb = (rand - 0.5)/(rand+0.1);
    ve = (rand - 0.5)/(rand+0.1);
    pb = (rand - 0.5)/(rand+0.1);
    pe = (rand - 0.5)/(rand+0.1);
    
    T1 = (pb-pe)/(0.5*vb+0.5*ve);
    T2 = (pb-pe)/(0.6*vb+0.4*ve);
    T3 = (pb-pe)/(0.4*vb+0.6*ve);

    if((T2-T1)*(T3-T1) > 0 && sign(T1) * sign(T2) > 0 && sign(T1) * sign(T3) > 0)
        error('sss1')
    elseif((sign(T1) * sign(T2) < 0 || sign(T1) * sign(T3) < 0) && (T1-T2)*(T1-T3) < 0)
        error('sss2')
    end

end

%%
clear
clc
vmax = 2.8;
amax = 3.2;
pe = (rand(1) - 0.5)/(rand(1)+0.01);
pb = 0;
for i=1:1e4
    vb = (rand*2-1)*vmax;
    ve = (rand*2-1)*vmax;
    
    T1 = (2*pe-2*pb)/(vb+ve);
    T2 = (3*pe-3*pb)/(2*vb+ve);
    T3 = (3*pe-3*pb)/(vb+2*ve);
    
    A1 = (- vb^2 - vb*ve - 3*vmax*vb - ve^2 - 3*vmax*ve);
    B1 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve - 6*pb*vmax + 6*pe*vmax);
    C1 = (- 9*pb^2 + 18*pb*pe - 9*pe^2);

    A2 = (- vb^2 - vb*ve + 3*vmax*vb - ve^2 + 3*vmax*ve);
    B2 = (6*pe*vb - 6*pb*ve - 6*pb*vb + 6*pe*ve + 6*pb*vmax - 6*pe*vmax);
    C2 = (- 9*pb^2 + 18*pb*pe - 9*pe^2);

    if(B1*B1-4*A1*C1 > 0)
        T_candidate1 = (-B1-sqrt(B1*B1-4*A1*C1))/(2*A1);
        T_candidate2 = (-B1+sqrt(B1*B1-4*A1*C1))/(2*A1);
        T4 = min([T_candidate2,T_candidate1]);
        T5 = max([T_candidate2,T_candidate1]);
    else
        T4 = -1;
        T5 = -1;
    end
    
    if(B2*B2-4*A2*C2 > 0)
        T_candidate1 = (-B2-sqrt(B2*B2-4*A2*C2))/(2*A2);
        T_candidate2 = (-B2+sqrt(B2*B2-4*A2*C2))/(2*A2);
        T6 = min([T_candidate2,T_candidate1]);
        T7 = max([T_candidate2,T_candidate1]);
    else
        T6 = -1;
        T7 = -1;
    end

%     D1 = -A1*(vb+ve);
%     D2 = A2*(vb+ve);
%     
%     if(D1 > 0 && D2 > 0)
%         Tmin_available = max([0.0, T1,T4,T5,T6,T7]);
%     elseif(D1 < 0 && D2 > 0)
%         Tk_a = [0,T1,T4,T5];
%         Tk_b = [0,T1,T6,T7];
%         Tk_a = sort(Tk_a);
%         Tk_b = sort(Tk_b);
% 
%         if(Tk_a(4) > Tk_b(4) && Tk_a(4) > max(T2,T3))
%             Tmin_available = max([0.0, Tk_a(3), Tk_b(4), T2, T3]);
%         else
%             Tmin_available = max([0.0, T2, T3]);
%         end
%     else
%         Tk_a = [0,T1,T4,T5];
%         Tk_b = [0,T1,T6,T7];
%         Tk_a = sort(Tk_a);
%         Tk_b = sort(Tk_b);
% 
%         if(Tk_b(4) > Tk_a(4) && Tk_b(4) > max(T2,T3))
%             Tmin_available = max([0.0, Tk_b(3), Tk_a(4), T2, T3]);
%         else
%             Tmin_available = max([0.0, T2, T3]);
%         end
%     end

    T_min_candidate = sort([T4,T5,T6,T7]);
    
    for i=1:4
        T = T_min_candidate(5-i);
        m = (T*(3*pb - 3*pe + 2*T*vb + T*ve))/(3*(2*pb - 2*pe + T*vb + T*ve));
        if((m < T && m > 0) || T < 0.0)
            break;
        else
            T = 0;
        end
    end

    % - amax*T^2 + (2*vb + 4*ve)*T + 6*pb - 6*pe == 0
% amax*T^2 + (2*vb + 4*ve)*T + 6*pb - 6*pe == 0
% - amax*T^2 + (- 4*vb - 2*ve)*T - 6*pb + 6*pe == 0
% amax*T^2 + (- 4*vb - 2*ve)*T - 6*pb + 6*pe == 0
    coes = [
        -amax, (2*vb + 4*ve), 6*pb - 6*pe;
        amax, (2*vb + 4*ve), 6*pb - 6*pe;
        -amax, -(4*vb + 2*ve), -6*pb + 6*pe;
        amax, -(4*vb + 2*ve), -6*pb + 6*pe;
        ];
    for i =1:4
        A = coes(i,1);
        B = coes(i,2);
        C = coes(i,3);

        if(sqrt(B*B-4*A*C) > 0)
            T = max(T, (-B+sqrt(B*B-4*A*C))/(2*A));
            T = max(T, (-B-sqrt(B*B-4*A*C))/(2*A));
        end
    end
    
%     T
    % 至此已经计算完毕，下面验证T是否满足要求

%     T = Tmin_available;
%     T = T_min_candidate(4);

%     m = (T*(3*pb - 3*pe + 2*T*vb + T*ve))/(3*(2*pb - 2*pe + T*vb + T*ve));

%     if(m > T || m < 0)
%         error('failed')
%     end
    
    M = [0,0,0,1
     0,0,1,0
     T^3, T^2, T,1
     3*T^2, 2*T,1,0];

    b = [pb;vb;pe;ve];

    x = M\b;

    A = x(1);
    B = x(2);
    C = x(3);
    D = x(4);

    k2 = 3*A;
    k1 = 2*B;
    k0 = C;

    m  = -k1/2/k2;
    
    vm = 3*A*m^2 + 2*B*m + C;
    a0 = 2*B;
    aT = 6*A*T+2*B;
    [vm,a0,aT]

%     if(m < T && m > 0)
%         vm = 3*A*m^2 + 2*B*m + C;
%         if((vm > vmax+1e-10 && vm < -vmax - 1e-10 || (abs(abs(vm) - vmax) > 1e-10)) ...
%                 &&)
%             error('sss');
%         end
%         
%     end

%     if(T1 <0  && T4 < 0 && T6 < 0 )
%         error('sss1')
%     end
% && ((D1 > 0 && D2 < 0) || (D1 < 0 && D2 > 0))

%     if(D1 > 0 && D2 < 0)
%         error('sss1')
%     end

%     if(D1 > 0 && D2 > 0)
%         error('sss2')
%     end
% 
%     if(D1 < 0 && D2 < 0)
%         error('sss3')
%     end

%     if(D1 > 0 && D2 > 0)
%         error('sss4')
%     end

end
%%

% T = T_min_candidate(2)

M = [0,0,0,1
     0,0,1,0
     T^3, T^2, T,1
     3*T^2, 2*T,1,0];

b = [pb;vb;pe;ve];

x = M\b;

A = x(1);
B = x(2);
C = x(3);
D = x(4);

k2 = 3*A;
k1 = 2*B;
k0 = C;


m  = -k1/2/k2

k2 * m^2 + k1 * m + k0


t = 0:0.01:T;
plot(t,k2*t.^2 + k1*t+k0)
