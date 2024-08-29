%%
%
% 下推 从 a0 v0 p0 -> ae ve pe 的s曲线计算
% 在 vmax amax jmax的约束下
%
% s曲线最多有 7 段：
%              j         a               v       
% 1 加加速段   jmax    a0 -> a_max     v0 -> vmax
% 2 匀加速段   0         a_max         v0 -> vmax  
% 3 减加速段  -jmax  amax -> 0         v0 -> vmax
% 4 匀速段     0          0             vmax 
% 5 加加速段   jmax     0 -> a_max   vmax -> ve
% 6 匀加速段   0         a_max       vmax -> ve
% 7 减加速段  -jmax  amax -> ae      vmax -> ve
%
% 每段运行时间为 T1 ... T7 每段结束后的状态为 ai vi pi，每段的jerk 为 ji
% 在理想情况下，对 j 积分：
% a1 = a0 + j1 * T1
% v1 = 
% 
% 可能有3种情况：
% CASE 1 ：可以达到 vmax 并匀速运动一段时间
% CASE 2 ：可以达到 -vmax 并匀速运动一段时间
% CASE 3 ：匀速段时间为0，且匀速段速度 为某个值 v4
%
% 在 CASE 1 下：
% 有：
% T1 = (amax - a0) / jmax
% T2 = (vmax - v0) / amax 
% T3 = amax / jmax


%%
clear
syms a0 v0 p0 ae ve pe T1 T2 T3 T4 T5 T6 T7 t jmax amax vmax



j1 = jmax
j2 = 0
j3 = -jmax
j4 = 0
j5 = -jmax
j6 = 0
j7 = jmax

a1(t) = int(j1,t) + a0
a2(t) = amax
a3(t) = int(j3,t) + a2(T2)
a4(t) = sym(0)
a5(t) = int(j5,t) + a4(T4)
a6(t) = -amax
a7(t) = int(j7,t) + a6(T6)
ae_e  = a7(T7)

v1(t) = int(a1,t) + v0
v2(t) = int(a2,t) + v1(T1)
v3(t) = int(a3,t) + v2(T2)
v4(t) = vmax
v5(t) = int(a5,t) + v4(T4)
v6(t) = int(a6,t) + v5(T5)
v7(t) = int(a7,t) + v6(T6)
ve_e  = v7(T7)

p1(t) = int(v1,t) + p0
p2(t) = int(v2,t) + p1(T1)
p3(t) = int(v3,t) + p2(T2)
p4(t) = int(v4,t) + p3(T3)
p5(t) = int(v5,t) + p4(T4)
p6(t) = int(v6,t) + p5(T5)
p7(t) = int(v7,t) + p6(T6)
pe_e  = p7(T7)


%%
% CASE 1：可以加速完整走出7段，即可达 amax vmax -amax
% 

% 因为可以达到 amax
T1 = (-a0 + amax) / jmax;
T3 = amax / jmax;
T2 = -((jmax*T1^2)/2 + a0*T1 + v0 - vmax + T3*amax - (T3^2*jmax)/2)/amax;


% 因为可以达到 -amax
T5 = amax / jmax;
T7 = (ae + amax) / jmax;
T6 = -((jmax*T5^2)/2 + ve - vmax + T7*amax - (T7^2*jmax)/2)/amax;


% 因为可以达到 vmax
T4 = -(p0 - pe + T6*(vmax - (T5^2*jmax)/2) + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + T5*vmax + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 - (T6^2*amax)/2 - (T7^2*amax)/2 - T7*((jmax*T5^2)/2 - vmax + T6*amax) + (T1^3*jmax)/6 - (T3^3*jmax)/6 - (T5^3*jmax)/6 + (T7^3*jmax)/6)/vmax;


%%
clear
p0 = 0.1;
v0 = -0.3;
a0 = 0.5;

pe = 100;
ve = 1.0;
ae = 0.2;

vmax = 2;
amax = 5;
jmax = 50;

% 因为可以达到 amax
T1 = (-a0 + amax) / jmax;
T3 = amax / jmax;
T2 = -((jmax*T1^2)/2 + a0*T1 + v0 - vmax + T3*amax - (T3^2*jmax)/2)/amax;


% 因为可以达到 -amax
T5 = amax / jmax;
T7 = (ae + amax) / jmax;
T6 = -((jmax*T5^2)/2 + ve - vmax + T7*amax - (T7^2*jmax)/2)/amax;


% 因为可以达到 vmax
T4 = -(p0 - pe + T6*(vmax - (T5^2*jmax)/2) + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + T5*vmax + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 - (T6^2*amax)/2 - (T7^2*amax)/2 - T7*((jmax*T5^2)/2 - vmax + T6*amax) + (T1^3*jmax)/6 - (T3^3*jmax)/6 - (T5^3*jmax)/6 + (T7^3*jmax)/6)/vmax;

T = T1 + T2 + T3 + T4 + T5 + T6 + T7;

dt = 0.01;
t_range = 0:dt:T;
p = zeros(size(t_range));

for i = 1:length(t_range)
    if(t_range(i) < T1)
        t = t_range(i);
        p(i) = (jmax*t^3)/6 + (a0*t^2)/2 + v0*t + p0;
    elseif(t_range(i) < T1 + T2)
        t = t_range(i) - T1;
        p(i) = p0 + T1*v0 + t*((jmax*T1^2)/2 + a0*T1 + v0) + (T1^2*a0)/2 + (T1^3*jmax)/6 + (amax*t^2)/2;
    elseif(t_range(i) < T1 + T2 + T3)
        t = t_range(i) - T1 - T2;
        p(i) = p0 + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + t*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + (T1^2*a0)/2 + (T2^2*amax)/2 + (T1^3*jmax)/6 + (amax*t^2)/2 - (jmax*t^3)/6;
    elseif(t_range(i) < T1 + T2 + T3 + T4)
        t = t_range(i) - T1 - T2 - T3;
        p(i) = p0 + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + t*vmax + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 + (T1^3*jmax)/6 - (T3^3*jmax)/6;
    elseif(t_range(i) < T1 + T2 + T3 + T4 + T5)
        t = t_range(i) - T1 - T2 - T3 - T4;
        p(i) = p0 + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + T4*vmax + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + t*vmax + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 + (T1^3*jmax)/6 - (T3^3*jmax)/6 - (jmax*t^3)/6;
    elseif(t_range(i) < T1 + T2 + T3 + T4 + T5 + T6)
        t = t_range(i) - T1 - T2 - T3 - T4 - T5;
        p(i) = p0 + T2*((jmax*T1^2)/2 + a0*T1 + v0) + t*(vmax - (T5^2*jmax)/2) + T1*v0 + T4*vmax + T5*vmax + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 + (T1^3*jmax)/6 - (T3^3*jmax)/6 - (T5^3*jmax)/6 - (amax*t^2)/2;
    elseif(t_range(i) < T1 + T2 + T3 + T4 + T5 + T6 + T7)
        t = t_range(i) - T1 - T2 - T3 - T4 - T5 - T6;
        p(i) = p0 + T6*(vmax - (T5^2*jmax)/2) + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + T4*vmax + T5*vmax + T3*((jmax*T1^2)/2 + a0*T1 + v0 + T2*amax) + (T1^2*a0)/2 + (T2^2*amax)/2 + (T3^2*amax)/2 - (T6^2*amax)/2 + (T1^3*jmax)/6 - (T3^3*jmax)/6 - (T5^3*jmax)/6 - (amax*t^2)/2 - t*((jmax*T5^2)/2 - vmax + T6*amax) + (jmax*t^3)/6;
    end

end


%%
%达不到 vmax 时

% clear


plot(diff(diff(p))/dt/dt)


%% 只有三段轨迹时
clear
syms a0 v0 p0 ae ve pe T1 T2 T3 jmax ac t



j1 = jmax
j2 = 0
j3 = -jmax

a1(t) = int(j1,t) + a0
a2(t) = ac
a3(t) = int(j3,t) + a2(T2)
ae_e  = a3(T3)

v1(t) = int(a1,t) + v0
v2(t) = int(a2,t) + v1(T1)
v3(t) = int(a3,t) + v2(T2)
ve_e(T1,T2,T3) = v3(T3)

p1(t) = int(v1,t) + p0
p2(t) = int(v2,t) + p1(T1)
p3(t) = int(v3,t) + p2(T2)
pe_e(T1,T2,T3)  = p3(T3)


% solve(a3(T3) == ae, T3)
T3 = (ac - ae)/jmax
% solve(a1(T1) == ae, T1)
T1 = -(a0 - ac)/jmax
% collect(solve(ve_e(T1,T2,T3) == ve, T2), ac)
T2 = (a0^2 - 2*ac^2 + ae^2 - 2*jmax*v0 + 2*jmax*ve)/(2*jmax*ac)

% collect(pe_e(T1,T2,T3)-pe, ac)    gives:
A = (- 6*a0^2 + 6*ae^2 + 12*jmax*v0 + 12*jmax*ve);
B = (24*jmax^2*p0 - 24*jmax^2*pe + 8*a0^3 - 8*ae^3 - 24*a0*jmax*v0 - 24*ae*jmax*ve);
C = - 3*a0^4 + 12*a0^2*jmax*v0 + 3*ae^4 + 12*ae^2*jmax*ve - 12*jmax^2*v0^2 + 12*jmax^2*ve^2;

A*ac^2 + B * ac + C == 0




%% test
clear
p0 = -0.1;
v0 = 0.3;
a0 = -0.2;

pe = -1;
ve = -1.0;
ae = -0.3;

jmax = -1;


A = (- 6*a0^2 + 6*ae^2 + 12*jmax*v0 + 12*jmax*ve);
B = (24*jmax^2*p0 - 24*jmax^2*pe + 8*a0^3 - 8*ae^3 - 24*a0*jmax*v0 - 24*ae*jmax*ve);
C = - 3*a0^4 + 12*a0^2*jmax*v0 + 3*ae^4 + 12*ae^2*jmax*ve - 12*jmax^2*v0^2 + 12*jmax^2*ve^2;


ac = (-B + sqrt(B*B-4*A*C))/(2*A)
% ac = (-B - sqrt(B*B-4*A*C))/(2*A)

T3 = (ac - ae)/jmax
T1 = -(a0 - ac)/jmax
T2 = (a0^2 - 2*ac^2 + ae^2 - 2*jmax*v0 + 2*jmax*ve)/(2*jmax*ac)

T= T1+T2+T3;

dt = 0.001;
t_range = 0:dt:T;
p = zeros(size(t_range));

for i = 1:length(t_range)
    if(t_range(i) < T1)
        t = t_range(i);
        p(i) = (jmax*t^3)/6 + (a0*t^2)/2 + v0*t + p0;
    elseif(t_range(i) < T1 + T2)
        t = t_range(i) - T1;
        p(i) = p0 + T1*v0 + t*((jmax*T1^2)/2 + a0*T1 + v0) + (T1^2*a0)/2 + (T1^3*jmax)/6 + (ac*t^2)/2;
    elseif(t_range(i) < T1 + T2 + T3)
        t = t_range(i) - T1 - T2;
        p(i) = p0 + T2*((jmax*T1^2)/2 + a0*T1 + v0) + T1*v0 + t*((jmax*T1^2)/2 + a0*T1 + v0 + T2*ac) + (T1^2*a0)/2 + (T2^2*ac)/2 + (T1^3*jmax)/6 + (ac*t^2)/2 - (jmax*t^3)/6;
    end

end

subplot(1,3,1)
plot(p)
subplot(1,3,2)
plot(diff(p)/dt)
subplot(1,3,3)
plot(diff(diff(p))/dt/dt)