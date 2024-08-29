%%
clear;
pt0 = 0;
vt0 = 0.5;
at0 = 0.0;

pa0 = 0;
va0 = 0;
aa0 = 0;

max_v = 1;
max_a = 2;
max_j = 5;

dt = 0.01;
%%
SIZE = 2000;
pt = zeros(1,2000);
vt = zeros(1,2000);
at = zeros(1,2000);

pa = zeros(1,2000);
va = zeros(1,2000);
aa = zeros(1,2000);

pt(1) = pt0;
vt(1) = vt0;
at(1) = at0;

pa(1) = pa0;
va(1) = va0;
aa(1) = aa0;



for i = 2:2000
    at(i) = at(i-1);
    vt(i) = vt(i-1) + at(i-1)*dt;
    pt(i) = pt(i-1) + vt(i-1)*dt + 0.5*at(i-1)*dt*dt;

    aa(i) = aa0;
    va(i) = va(i-1) + aa(i-1)*dt;
    pa(i) = pa(i-1) + va(i-1)*dt + 0.5*aa(i-1)*dt*dt;
end
%%
subplot(2,3,1)
plot(pt)
subplot(2,3,2)
plot(vt)
subplot(2,3,3)
plot(at)
subplot(2,3,4)
plot(pa)
subplot(2,3,5)
plot(va)
subplot(2,3,6)
plot(aa)

%%
clear
syms a b c d e f T t

func = a*t^5 + b*t^4 + c*t^3 + d*t^2 + e*t^1 + f

dfunc = diff(func, t)

%%
clear
syms T1 T2 T3 T4 T5 T6 T7
syms p0 v0 a0 pe ve ae T
syms max_v max_a max_j

a1 = a0 + max_j * T1;
v1 = v0 + 0.5*(a0+a1)*T1;
p1 = v0 * T1 + max_j * T1;











