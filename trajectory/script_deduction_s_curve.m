%%
syms v0 v1 T T1

% T1 = 2*sqrt((v0 - v1)/j)
eqn =  (v0+v1)/2 * T1 + v1 * (T-T1)
expand((v0+v1)/2 * T1 + v1 * (T-T1))

%% 
expand((v0-v1)^4)

%%
p0 = 0.1;
v0 = 1.0;
p1 = 2.5;

j = 0.8;
a = 0.5;

T = 1.9;

% s_curve(0,p0,v0,p1,v,a,j,T)

% v1 = 0.5;
% si = sign(v1 - v0);
% T1 = si/a*(v1 - v0) + a/j;


dt = 0.001;
t = 0:dt:T;
p = zeros(size(0:dt:T));

% s_curve(0,p0,v0,p1,a,j,T)

for i = 1:length(t)
%     if(t(i) < a/j)
%         p(i) = p0 + v0 * t(i) + si/6*j*t(i)^3;
%     elseif(t(i) < (T1 - a/j))
%         p(i) = p0 + v0 * a / j + si/6 * a^3 / j^2 + (v0 + si/2*a^2/j) * (t(i) - a/j)+ si*a/2*(t(i)-a/j)^2; 
%     elseif(t(i) < T1)
%         p(i) = p0 + (v0+v1)/2*T1-(v1*(T1-t(i)) - si/6*j*(T1-t(i))^3);
%     else
%         p(i) = p0 + (v0+v1)/2*T1 + v1 * (t(i) - T1);
%     end
    [p(i), v1] = s_curve(t(i),p0,v0,p1,a,j,T);
end

hold on
plot(t,p);
plot(t(2:end) - dt/2,diff(p)/dt);
plot(t(3:end) - dt/2*3,diff(diff(p(1:end)))/dt/dt);

axis equal

%%
p0 = 0.1;
v0 = 1.0;
p1 = 2.1;

v = 0.5;
j = 0.8;
a = 0.5;

T = 3;

% s_curve(0,p0,v0,p1,v,a,j,T)

% v1 = 0.5;
% T1 = 2 * sqrt( abs(v0 - v1) / j );
% p1 = p0 + (v0+v1)/2 * T1 + v1 * (T - T1);



dt = 0.001;
t = 0:dt:T;
p = zeros(size(0:dt:T)) ;

for i = 1:length(t)
%     if(t(i) < T1 / 2)
%         p(i) = p0 + v0 * t(i) + sign(v1 - v0) * 1/6*j*t(i)^3;
%     elseif(t(i) < T1)
%         p(i) = p0 + (v0+v1)/2*T1-(v1 * (T1-t(i)) - sign(v1 - v0) * 1/6*j*(T1-t(i))^3);
%     else
%         p(i) = p0 + (v0+v1)/2*T1 + v1 * (t(i) - T1);
%     end
    [p(i), v1] = s_curve(t(i),p0,v0,p1,a,j,T);
end

hold on
plot(t,p);
plot(t(2:end) - dt/2,diff(p)/dt);
plot(t(3:end) - dt/2*3,diff(diff(p))/dt/dt);

axis equal
%%
s_curve(0,p0,v0,p1,v1,a,j,T)

%%
% T*v1 + (T1*v0)/2 - (T1*v1)/2 = p
v1 = sign(v1-v0)*j/4*T1^2+v0
T*v1 + (T1*v0)/2 - (T1*v1)/2 - p1



%%
% sig is sign(v1 - v0)
clear
syms sig p1 T T1 v0 j a
v1 = sig*j/4*T1^2+v0
collect((T*v1 + (T1*v0)/2 - (T1*v1)/2 - p1)*10, T1)



%%
syms p0 v0 si T1 a j T
v1 = si*a*T1 + (v0 - si*a^2/j)
% p0 + (v0+v1)/2 * T1 + v1 * (T-T1)
collect(p0 + (v0+v1)/2 * T1 + v1 * (T-T1), T1)
%% 
syms sig p1 T v1 v0
expand((p1 - T*v1)^2 - sig/j*(v0-v1)^2)
collect((p1 - T*v1)^2 - sig/j*(v0-v1)^2, v1)



%% 
syms a b x
root(x^3 - a*x^2 + b)

solve(x^3 - a*x^2 + b == 0, x)

%%
a=rand(1,1)-0.5;
b=rand(1,1)-0.5;
c=rand(1,1)-0.5;
d=rand(1,1)-0.5;

clc
x = cubic_equation_solve(a,b,c,d)
a*x.^3 + b*x.^2 + c * x + d
%%
syms va a j vb
t1 = (va-vb)/a+a/j;
p1 = t1*(va+vb)/2;
collect(-p1*2*a, vb)

t2 = 2 * sqrt( (va-vb) / j );
p2 = t2*(va+vb)/2;
collect(p2^2, vb)
%%
syms va a j vb
vb = - a^2/j + va
t1 = (va-vb)/a+a/j;
p1 = t1*(va+vb)/2
expand(p1)
%%
syms v1 a j v2 v

v = v1+a^2/j
t1 = (v-v1) / a + a/j
p1 = (v+v1)/2*t1
t2 = (v-v2) / a + a/j
p2 = (v+v2)/2*t2

T = t1+t2
% p = p1+p2
% collect(T,v)
expand(T)

%%
syms v1 a j v2 v T

% v = v2+a^2/j
t1 = 2 * sqrt((v-v1) / j )
t2 = 2 * sqrt((v-v2) / j )

eq1=(t1+t2)^2-T^2

eq1 = (2*((v - v1)/j)^(1/2) + 2*((v - v2)/j)^(1/2))^2 - T^2

term1 = 8*(v/j - v1/j)^(1/2)*(v/j - v2/j)^(1/2) 
term2 = (8*v)/j - (4*v1)/j - (4*v2)/j - T^2
term1^2 - term2^2

eq2 = 64*(v/j - v1/j)*(v/j - v2/j) - (T^2 - (8*v)/j + (4*v1)/j + (4*v2)/j)^2
expand(eq2*j^2)
collect(eq2*j^2,v)

% expand(a^2*((T - a/j + v2/a)^2 + (4*v1)/j))
%%
syms va vb vb_max Ta Tb T pt a j
Tb = (vb_max-vb) / a + a/j

eq = Ta * (va + vb_max)/2 + vb_max * (T-Ta-Tb) + Tb*(vb_max + vb)/2 - pt
collect(eq,vb)

%%
clear
syms max_v a j Tb T Ta pa pt
% Tb = 2*a/j
% vb = max_v - j*Tb*Tb/4
vb = max_v - Tb*a + a^2/j
pb = (vb + max_v)/2*Tb
pmid = (T-Ta-Tb)*max_v

eq = pa + pb + pmid - pt

collect(-eq*2/a, Tb)
%%
clear
syms v a j Tb T Ta pacc_max pt
% Tb = 2*a/j
% vb = max_v - j*Tb*Tb/4
vb = v - Tb*a + a^2/j
pb = (vb + v)/2*Tb

eq = pacc_max + Tb*(a^2/(2*j) - (Tb*a)/2 + v) - pt

collect(eq,Tb)

%% 末端速度为0，a段达不到最大加速度，b段达不到最大加速度
clear
syms va j Ta pt T f(Ta)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = sqrt(4*va/j + Ta*Ta)
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l  = la+lb+lc
% collect(l,Ta)

func_la_lb_lc_pt = ...
((- (j*Ta^2)/8 - va/2)*Tb)^2 - ... 
    (pt - Ta*((j*Ta^2)/8 + va) ...
    - ((j*Ta^2)/4 + va)*(T - Ta))^2;
% collect(func_la_lb_lc_pt, Ta);
% 
f(Ta) = func_la_lb_lc_pt;
((j*Ta^2)/8 + va/2)^2*(Ta^2 + (4*va)/j) - (Ta*((j*Ta^2)/8 + va) - pt + ((j*Ta^2)/4 + va)*(T - Ta))^2


diff(lb,Ta)

K  = (j*Ta^2) + 4*va

eq = (Ta*K) / (8*Tb) +(Ta*j*Tb)/4 - diff(lb,Ta)
eq1 = expand(eq)

eq = - (Ta/Tb + 1)*(K/4) - (Ta*j*(Ta - T + Tb))/2 ...
    - diff(lc,Ta)
eq1 = expand(eq)

clear
syms j Ta va K Tb T
dl = (3*j*Ta^2)/8 + va + (Ta*K) / (8*Tb) +(Ta*j*Tb)/4 ...
        - (Ta/Tb + 1)*(K/4) - (Ta*j*(Ta - T + Tb))/2

expand(dl)
%% 末端速度为0，a段达不到最大加速度，b段可以达到最大加速度
clear
syms va j Ta pt T a f(Ta)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = (v + a^2/j)/a
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l  = la+lb+lc-pt
% collect(l,Ta)
eq = va - v*((Ta*j)/2 + 1) + (3*Ta^2*j)/8 + (Ta*j*(a^2/j + v))/4 - (Ta*j*(Ta + a^2/j - T + v))/2 + (Ta*j*v)/4
eq1 = eq - diff(l,Ta);
expand(eq1);
expand(eq)
% func_la_lb_lc_pt = ...
% ((- (j*Ta^2)/8 - va/2)*Tb)^2 - ... 
%     (pt - Ta*((j*Ta^2)/8 + va) ...
%     - ((j*Ta^2)/4 + va)*(T - Ta))^2;
% collect(func_la_lb_lc_pt, Ta);
% 
% f(Ta) = func_la_lb_lc_pt;
%% 末端速度为0，a段可以达到最大加速度，b段可以达到最大加速度
clear
syms va j Ta pt T a f(Ta)
v  = va + Ta*a - a^2/j
la = (v+va)/2*Ta
Tb = (v + a^2/j)/a
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l  = la+lb+lc-pt
f(Ta) = l
% collect(l,Ta)
eq = va - v*((Ta*j)/2 + 1) + (3*Ta^2*j)/8 + (Ta*j*(a^2/j + v))/4 - (Ta*j*(Ta + a^2/j - T + v))/2 + (Ta*j*v)/4
eq1 = eq - diff(l,Ta);
expand(eq1);
expand(eq)
% func_la_lb_lc_pt = ...
% ((- (j*Ta^2)/8 - va/2)*Tb)^2 - ... 
%     (pt - Ta*((j*Ta^2)/8 + va) ...
%     - ((j*Ta^2)/4 + va)*(T - Ta))^2;
% collect(func_la_lb_lc_pt, Ta);
% 
% f(Ta) = func_la_lb_lc_pt;

%%
clear;
syms T v a j 
T=2*a/j
j*T*T/4

% collect((pt-la)^2 - lb^2, Ta)
%% 末端速度为0
clear
syms v va j Ta pt T a
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = v/ a + a/j;
vb = 0
lb = Tb*(v+vb)/2


collect(la + lb - pt, Ta)
%% subcase 1
clear
syms v va j Ta pt T a
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = T - Ta;
vb = v - j*Tb*Tb/4
lb = Tb*(v+vb)/2


collect(-la - lb + pt, Ta)
%% subcase 2
clear
syms v va j Ta pt T a
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = T - Ta;
vb = v - Tb*a + a^2/j
lb = Tb*(v+vb)/2


collect(la + lb - pt, Ta)
%% subcase 4
clear
syms v va j Ta pt T a
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = T - Ta;
vb = v - Tb*a + a^2/j
lb = Tb*(v+vb)/2


collect(la + lb - pt, Ta)
%% subcase 3
syms v va j Ta pt T a
v  = va + Ta*a - a^2/j
la = (v+va)/2*Ta
Tb = T - Ta;
vb = v - j*Tb*Tb/4
lb = Tb*(v+vb)/2

collect(la + lb - pt, Ta)


%% 有匀速段，可降速为0
clear
va = 0.1;
a = 5;
j = 5;



Ta = 0.1
v  = s_acc_vend(va,a,j,Ta)
Tb = s_acc_time(0,v,a,j)
Tc = 0.2
T = Ta+Tb+Tc

la = Ta*(va+v)/2;
lb = Tb*v/2;
lc = Tc*v;

pt = la+lb+lc;
% v  = va + j*Ta*Ta/4
% la = (v+va)/2*Ta
% Tb = sqrt(4*va/j + Ta*Ta)
% lb = Tb*v/2
% 
% la = (v+va)/2*Ta
% lb = Tb*v/2


k0 = va^3/j - (pt - T*va)^2;
k2 = (3*va^2)/4 + (T*j*(pt - T*va))/2;
k3 = (j*(T*va-pt))/4;
k4 = (3*j*va)/16 - (T^2*j^2)/16;
k5 = (T*j^2)/16;
k5*Ta^5 +k4*Ta^4 +k3*Ta^3 + k2*Ta^2 + k0
%% 无匀速段
clear
va = 0.1;
vb = 0.05;
a = 5;
j = 5;


Ta = 0.1
v  = s_acc_vend(va,a,j,Ta)
Tb = s_acc_time(vb,v,a,j)
Tc = 0.0
T = Ta+Tb+Tc

la = Ta*(va+v)/2;
lb = Tb*(vb+v)/2;

pt = la+lb

k2 = T*j/8
k1 = -3*T^2*j/8
k0 = pt - (T*(2*va - (T^2*j)/4))/2

k2*Ta*Ta + k1*Ta + k0

roots([k2,k1,k0])
%% 有匀速段
clear
va = 1.2;
vb = 0.05;
a = 5;
j = 5;


Ta = 1.98
v  = s_acc_vend(va,a,j,Ta)
Tb = s_acc_time(vb,v,a,j)
Tc = 0.0
T = Ta+Tb+Tc

la = Ta*(va+v)/2;
lb = Tb*(vb+v)/2;

pt = la+lb

k3 = - j/8
k2 = (T*j)/4 - a/2
k1 = - a^2/(2*j) + T*a
k0 = (T*(a^2/j - T*a + 2*va))/2 - pt

k3*Ta*Ta*Ta + k2*Ta*Ta + k1*Ta + k0

roots([k3,k2,k1,k0])


s_make_s_curve(0,va,pt,100,a,j,T)




%%
syms v va j Ta pt T
Tb = 2*sqrt(va/j + Ta*Ta)
collect(-(T-Ta)^2+Tb^2,Ta)

A=3
B=2*T
C=(4*va)/j-T^2

Ta = (-B+sqrt(B^2-4*A*C))/4*A
Tb = T - Ta
collect(-(T-Ta)^2+Tb^2,Tb)
collect(expand(-(T-Ta)^2+Tb^2),T)

la = j/8*Ta^3 + va*Ta + Tb*v/2

% (-B-sqrt(B^2-4*A*C))/4*A
%% 
clear
syms v va j Ta pt T a
A=j/4/a
B=1
C=- T + a/j + va/a

Ta = (-B+sqrt(B^2-4*A*C))/4*A
2*sqrt(a^2/j^2 - va/j)

%% 
clear
syms a v j Ta Tb v1 v2 T
t1 = 2 * sqrt((v-v1) / j )
t2 = 2 * sqrt((v-v2) / j )

expand((t1+t2)^2 -T^2)
A = 8*(v/j - v1/j)^(1/2)*(v/j - v2/j)^(1/2)
B = (8*v)/j - (4*v1)/j - (4*v2)/j - T^2

collect(A^2-B^2,v)


((16*(T^2 + (4*v1)/j + (4*v2)/j))/j - (64*v1)/j^2 - (64*v2)/j^2)*v + (64*v1*v2)/j^2 - (T^2 + (4*v1)/j + (4*v2)/j)^2
solve(A^2-B^2,v)

%%
clear
syms va vb a j pt
T = (vb-va) / a + a/j
l = T * (va+vb)/2
coeffs(l,vb)

k2 = 1/(2*a);
    k1 = a/(2*j);
    k0 = (va*(a/j - va/a))/2;

syms f(vb)
f(vb) = l
expand(f(vb))
collect(f(vb),vb)
%%
syms v1 v2 a j pt v;
T1 = 2 * sqrt((v-v1) / j )
T2 = 2 * sqrt((v-v2) / j )
l = T1 * (v1 + v) / 2 + T2 * (v2 + v) / 2

diff(l,v) -(...
    T1/2 + T2/2 + (v + v1)/(j*T1) + (v + v2)/(j*T2)... 
)
% eq = expand(l^2 - pt^2)
% collect(eq,v)
% coeffs(eq,v)

% deq = diff(eq,v)
% collect(deq,v)
% solve(deq,v)
% 
% f(v) = eq
% 
% value = f(v1/2 - v2/6)
% 
% expand(value)
%%
syms v1 v2 a j pt v f(v);
T1 = 2 * sqrt((v-v1) / j )
T2 = (v-v2) / a + a/j
l1 = T1 * (v1 + v) / 2
l2 = T2 * (v2 + v) / 2


l = T1 * (v1 + v) / 2 + T2 * (v2 + v) / 2
eq = expand(l1^2 - (pt-l2)^2)
collect(eq,v)
coeffs(eq,v)

eq = l - pt
deq = diff(eq,v)
collect(deq,v)

deq - T1/2 - T2/2 - (v + v2)/(2*a) - (v+v1)/(j*T1)
%%
syms v1 v2 a j pt v;
T1 = (v-v1) / a + a/j
T2 = (v-v2) / a + a/j
l1 = T1 * (v1 + v) / 2
l2 = T2 * (v2 + v) / 2
l  = l1 + l2

f(v) = l-pt;

df = diff(f,v)
% collect(l-pt,v)
% coeffs(l-pt,v)
% solve(l-pt,v)
%%
syms k3 k2 k1 k0 x z
solve(x^4+k3*x^3 + k2 * x^2 + k1 * x + k0,x)
%%
clear
syms va j Ta pt T f(Ta) g(T,Ta,va,j)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = sqrt(4*va/j + Ta*Ta)
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l  = la+lb+lc

term = ((j*Ta^2)/4 + va)

rhs = (pt - Ta*((j*Ta^2)/8 + va)) / ((j*Ta^2)/4 + va)
lhs = (l - Ta*((j*Ta^2)/8 + va))  / ((j*Ta^2)/4 + va)

g(T,Ta,va,j) = lhs
g(0.4,0.3,2.0,0.2)

syms t
term = t
lhs1 = -(term*(Ta - T + (Ta^2 + (4*va)/j)^(1/2)) - (term*(Ta^2 + (4*va)/j)^(1/2))/2)/term

lhs2 = (2*T - 2*Ta - ((j*Ta^2 + 4*va)/j)^(1/2))/2


g(T,Ta,va,j) = lhs2
result = g(0.4,0.3,2.0,0.2)

% f(0.4,0.3,2.0,0.8)


% va = 0.8
% j  = 2.0
% Ta = 0.2
% T  = 0.4

%% 
clear
% syms va j Ta pt T f(Ta) g(T,Ta,va,j)
% v  = va + j*Ta*Ta/4
% la = (v+va)/2*Ta
% Tb = sqrt(4*va/j + Ta*Ta)
% lb = Tb*v/2
% lc = (T - Ta - Tb)*v

syms T va j v pt f(v) g(T,Ta,va,j) lhs_f(T,v,va,j,pt) rhs(T,v,va,j,pt) 
Ta = sqrt((v-va)*4/j)
la = (v+va)/2*Ta
Tb = sqrt(4*va/j + Ta*Ta)
lb = Tb*v/2
lc = (T - Ta - Tb)*v

f(v) = la + lb + lc
g(T,va,j,v) = la + lb + lc

pt_v = g(2.8 , 1.2 , 0.8 , 2.0)

expand(la+lb+lc)

lhs1 = (v-va)*sqrt(v/j - va/j) + v*sqrt(v/j)
rhs1 = T*v - pt

lhs_f(T,v,va,j,pt) = lhs1;
lhs_v = lhs_f(2.8 , 1.2 , 0.8 , 2.0 , pt_v)
rhs_f(T,v,va,j,pt) = lhs1;
rhs_v = rhs_f(2.8 , 1.2 , 0.8 , 2.0 , pt_v)
lhs_v - rhs_v

lhs2 = (v-va)*sqrt(v-va) + v*sqrt(v)
rhs2 = (T*v - pt)*sqrt(j)

lhs_f(T,v,va,j,pt) = lhs2;
lhs_v = lhs_f(2.8 , 1.2 , 0.8 , 2.0 , pt_v)
rhs_f(T,v,va,j,pt) = lhs2;
rhs_v = rhs_f(2.8 , 1.2 , 0.8 , 2.0 , pt_v)
lhs_v - rhs_v
%%
syms b c d
q = (3.0*c - (b*b))/9.0;
r = (-(27.0*d) + b*(9.0*c - 2.0*(b*b)))/54;

disc = q*q*q + r*r;
term1 = (b/3.0);

s = r + sqrt(disc);
s = sign(s)*abs(s)^(1/3);
t = r - sqrt(disc);
t = sign(t)*abs(t)^(1/3);

x = -term1 + s + t;
%%
n = 1e5;

sa = sign(rand(1,n) - 0.5);
sb = sign(rand(1,n) - 0.5);
sc = sign(rand(1,n) - 0.5);
sd = sign(rand(1,n) - 0.5);

a = sa.*rand(1,n)*100;
b = sb.*rand(1,n)*100;
c = sc.*rand(1,n)*100;
d = sd.*rand(1,n)*100;

tic
for i =1:n
%     i
    x = cubic_equation_solve(a(i),b(i),c(i),d(i));

%     for j=1:length(x)
%         v = abs(a(i)*x(j)^3 + b(i)*x(j)^2+c(i)*x(j)+d(i));
%         r = roots([a(i),b(i),c(i),d(i)]);
%         
%         
% 
%         if(v > 1e-6 && min(abs(r - x(j))) > 1e-14 * abs(x(j)))
%             v
%             x
%             roots([a(i),b(i),c(i),d(i)])
%             min(abs(r - x(j)))
%             i
%             error("failed");
%         end
%     end
end
toc
%%
clc
x = cubic_equation_solve(a(i),b(i),c(i),d(i))
for j=1:length(x)
    v = abs(a(i)*x(j)^3 + b(i)*x(j)^2+c(i)*x(j)+d(i))
    if(v > 1e-10)
        v
        error("failed");

    end
end
%%
syms r x f(x) df(x) ddf(x) dddf(x) g(x) dg(x) ddg(x) dddg(x)
f(x) = (r + sqrt(r^2 + x^3))^(1/3)
df(x) = diff(f,x)
ddf(x) = diff(df)
dddf(x) = diff(ddf)

g(x) = (r + (r^2)^(1/2))^(1/3) + 1/((r + (r^2)^(1/2))^(2/3)*(r^2)^(1/2))*(1/6)*x.^3
dg(x) = diff(g,x)
ddg(x) = diff(dg,x)
dddg(x) = diff(ddg,x)

r=1;
g(0)
dg(0)
ddg(0)
dddg(0)
%%
syms r x f(x) df(x) ddf(x) dddf(x) g(x) dg(x) ddg(x) dddg(x)
f(x) = -(-r + sqrt(r^2 + x^3))^(1/3)
df(x) = diff(f,x)
ddf(x) = diff(df)
dddf(x) = diff(ddf)

h(x) = (r - (r^2)^(1/2))^(1/3) - 1/((r - (r^2)^(1/2))^(2/3)*(r^2)^(1/2))*(1/6)*x.^3
dh(x) = diff(h,x)
ddh(x) = diff(dh,x)
dddh(x) = diff(ddh,x)

r = 1
f(x) = -(-1 + sqrt(1 + x^3))^(1/3)
df(x) = diff(f,x)
ddf(x) = diff(df)
dddf(x) = diff(ddf)
vpa(df(1e-11))
% syms y
% f = (1-y^3)^2-1

%%
r = 1;
x=-0.1:1e-6:0.1;
plot(x,(r + sqrt(r^2 + x.^3)).^(1/3));
hold

g = (r + (r^2)^(1/2))^(1/3) + 1/((r + (r^2)^(1/2))^(2/3)*(r^2)^(1/2)) * (1/6)*x.^3;
plot(x,g)
h = 1.259921049894873 + 0.104993420824573*x.^3;
plot(x,h)

%%
r = 1;
x=0:1e-6:0.001;
t = -(-1 + sqrt(1 + x.^3)).^(1/3);
plot(x,-(-1 + sqrt(1 + x.^3)).^(1/3));
hold

% g = -0.79375086002259034701847785911451*x;
g = -0.79370052624990633724219768806982*x;
% g = -0.79370052598400404903731454005881*x;
% g = -0.79370070779435684530997705433948*x;

g(end) - t(end)
plot(x,g)



%%
    syms va vb a j pt
    T = 2 * sqrt((vb-va) / j )
    l = T * (va+vb)/2

%%
clear
syms va j Ta a T
v  = va + j*Ta*Ta/4
la = j/8*Ta^3 + va*Ta
Tb = T-Ta
vb = v - Tb*a + a^2/j;
lb = Tb*(v+vb)/2
l  = la + lb

%%
syms va j Ta a T
T = sqrt(va*8/3/j)
vb = va - j * T^2 / 4
%%
syms va j Ta pt T a f(Ta)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = v + a^2/j
lb = Tb*v/2
lc = (T - Ta - Tb)*v


la + lb + lc

%%
Ta_ = min([T, s_acc_time(va, max_v, a, j), 2*a/j])
v_  = va + j*Ta_*Ta_/4
Tb_ = T-Ta_
vb_ = max(v_ - j*Tb_*Tb_/4, 0)
l_  = Ta_ * (va + v_) / 2.0 + Tb_ * (vb_ + v_) / 2.0
%%
% Ta = min([T, s_acc_time(va, max_v, a, j)])
% v  = s_acc_vend(va,a,j,Ta)
% la = (va + v)*Ta/2
% Tb = T-Ta
% vb = s_acc_vend(v,-a,-j,Tb)
% lb = Tb*(v+vb)/2
% l  = la + lb
i = 1;

Tmax = zeros(1,n);
Tmin = zeros(1,n);
if(i==1)
    for j = 1:n
        [Tmax(j),Tmin(j)] = s_s_curve_Tmax_Tmin(0,0,pos(i,j),max_vb(i,j),vel(i,j),acc(i,j),jerk(i,j));
    end
else
    for j = 1:n
        v0 = [0        0.0287473827196604                         0                         0                   0 ];
        [Tmax(j),Tmin(j)] = s_s_curve_Tmax_Tmin(pos(i-1,j),v0(j),pos(i,j),max_vb(i,j),vel(i,j),acc(i,j),jerk(i,j));
    end
end
Tmax
Tmin
%%
clear
syms va j Ta pt T a f(Ta)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = (v + a^2/j)/a
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l = la+lb+lc

collect(l-pt,Ta)

%%
syms va T a j pt Ta 
    k4 = -j^2/(32*a);
    k3 = -j/8;
    k2 = (T*j)/4 - a/8 - (j*va)/(4*a);
    k0 = T*va - pt - va^2/(2*a) - (a*va)/(2*j);

syms f(Ta) g(Ta)
f(Ta) = k4*Ta^4 + k3*Ta^3 + k2*Ta^2 + k0
g(Ta) = 4*k4*Ta^3 + 3*k3*Ta^2 + 2*k2*Ta
solve(g,Ta)

% B^2 - 4AC
expand((3*k3)^2 - 4*(4*k4)*2*k2)
expand(- a^2/8 + (T*j)/4 - (j*va)/4 + 9/64)

syms r1 r2
r1 = (3*a + (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j);
r2 = (3*a - (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j);

% T = va/a+a/j
% r2 = (3*a - (a^2 + (T*j*a)*16 - (j*va)*16)^(1/2))/(2*j)

% expand(r1)
% expand(r2)
% expand(f(r1))
% expand(f(r2))

eq = f(0)
T = va/a+a/j

T*va - pt - va^2/(2*a) - (a*va)/(2*j)
expand(T*va - pt - va^2/(2*a) - (a*va)/(2*j))
%%
clear
syms va j pt T v vb f(Ta)
Ta = 2 * sqrt((v-va) / j )
Tb = 2 * sqrt((v - vb) / j )
eq = (Ta+Tb)^2
left = 8*(v/j - va/j)^(1/2)*(v/j - vb/j)^(1/2)
right = T^2 + (4*va)/j + (4*vb)/j - (8*v)/j
expand(left^2 - right^2)
collect(left^2 - right^2,v)
expand(((16*(T^2 + (4*va)/j + (4*vb)/j))/j - (64*va)/j^2 - (64*vb)/j^2))
expand((64*va*vb)/j^2 - (T^2 + (4*va)/j + (4*vb)/j)^2)

% la = Ta*(va+v)/2
% lb = Tb*(v+vb)/2
% lc = (T-Ta-Tb)*v
% l  = la+lb
%%
clear
syms va j pt T v vb f(Ta)
Ta = 2 * sqrt((v-va)/j)
Tb = 2 * sqrt(v/j)
eq = (Ta+Tb)^2
expand(eq)
left = 8*(v/j - va/j)^(1/2)*(v/j)^(1/2)
right = T^2 + (4*va)/j - (8*v)/j
expand(left^2 - right^2)
collect(left^2 - right^2,v)
expand((T^2 + (4*va)/j)^2)
expand(((16*(T^2 + (4*va)/j))/j - (64*va)/j^2))
solve(left^2 - right^2,v)
%%
clear
syms va j v a Ta Tb T f(v)
% Ta = (v-va)/a + a/j
% Tb = v/a + a/j
% 
% eq = expand(Ta^2 - (T-Tb)^2)
% collect(eq,v)
% f(v) = eq
% f(0)

solve(va + j*Ta*Ta/4 - a^2/j,Ta)
%%
% s_s_curve( ...
%                 t(j), ...
%                 p0(r), ...
%                 v0(r), ...
%                 pos(i,r), ...
%                 vend(i,r), ...
%                 vel(i,r), ...
%                 acc(i,r), ...
%                 jerk(i,r), ...
%                 T(i), ...
%                 Ta(i,r), ...
%                 Tb(i,r));
k=5;
Ta2 = Ta(k,2);
Tb2=Tb(k,2);
va2 = vb(k-1,2);
vb2=vb(k,2);
v2= real_v(k,2);
a2 = acc(k,2);
j2 = jerk(k,2);
T2 =T(k);
pa2 = pos(k-1,2);
pb2 = pos(k,2);

p_out = zeros(size(0:0.001:T2));
i=1;
for(t = 0:0.001:T2)
p_out(i) = s_s_curve( ...
                t, ...
                pa2, ...
                va2, ...
                pb2, ...
                vb2, ...
                v2, ...
                a2, ...
                j2, ...
                T2, ...
                Ta2, ...
                Tb2);
i = i+1;
end

plot(0:0.001:T2,p_out)
% plot(0:0.001:T2-0.001,diff(p_out))
% s_s_curve( ...
%                 T2, ...
%                 pa2, ...
%                 va2, ...
%                 pb2, ...
%                 vb2, ...
%                 v2, ...
%                 a2, ...
%                 j2, ...
%                 T2, ...
%                 Ta2, ...
%                 Tb2) - pa2

s_make_s_curve(pa2,va2,pb2,v2,a2,j2,T2)
s_curve_length(va2,vb2,v2,a2,j2,T2)


%%
clear
syms va Ta j a T pt
eq = (va + j*Ta*Ta/8)*Ta...
        + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
        + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4)

(va + j*Ta*Ta/8)*Ta...
            + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
            + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4) - pt

expand(eq)

expand(T*va - (Ta^3*j)/8 - (va/2 + Ta^2*j/8)*(Ta^2 + (4*va)/j)^(1/2) + (T*Ta^2*j)/4 -eq)




va = 0;
vb = 0;
% vc = 
% T=1.1128874161511026;
T=500.50581926029304;
a=50;
j=50;
pt = 1.9687164710040281;
% Ta = newton_raphson_binary_search(@(Ta)(...
%             (va + j*Ta*Ta/8)*Ta...
%             + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
%             + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4) - pt)...
%             ,0,2.0,10*eps)

Ta = newton_raphson_binary_search(@(Ta)(...
        T*va - (Ta^3*j)/8 ...
        - (va/2 + Ta^2*j/8)*(Ta^2 + 4*va/j)^(1/2) ...
        + (T*Ta^2*j)/4 - pt)...
            ,0,2.0,10*eps)

syms va j Ta pt T g(Ta)
g(Ta,j,va,T,pt) = T*va - (Ta^3*j)/8 ...
        - (va/2 + Ta^2*j/8)*(Ta^2 + 4*va/j)^(1/2) ...
        + (T*Ta^2*j)/4 - pt;

g(2,50,0,500,2.0)
%%
clear
syms va j Ta pt T f(Ta)
v  = va + j*Ta*Ta/4
la = (v+va)/2*Ta
Tb = sqrt(4*va/j + Ta*Ta)
lb = Tb*v/2
lc = (T - Ta - Tb)*v

l = la+lb+lc - pt

expand(l)

f(Ta,j,va,T,pt) = l

f(100,50,0,500,2.0)

% Ta = newton_raphson_binary_search(@(Ta)(...
%             (va + j*Ta*Ta/8)*Ta...
%             + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
%             + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4) - pt)...
%             ,0,2.0,10*eps)
% 
% 
% Ta = newton_raphson_binary_search(@(Ta)(...
%             (va + j*Ta*Ta/8)*Ta...
%             + sqrt(4*va/j + Ta*Ta)*(va + j*Ta*Ta/4)/2 ...
%             + (T - Ta - sqrt(4*va/j + Ta*Ta))*(va + j*Ta*Ta/4) - pt)...
%             ,0,2.0 ...
%             ,10*eps)


%%
clear;
syms va j a Ta T pt

v  = va + j*Ta*Ta/4
la = j/8*Ta^3 + va*Ta
Tb = T-Ta
vb = v  - j*Tb*Tb/4
lb = Tb*(v+vb)/2

expand(la + lb - pt)
%%
syms v1 v2 a j pt v;
T1 = 2 * sqrt((v-v1) / j )
T2 = 2 * sqrt((v-v2) / j )
l = T1 * (v1 + v) / 2 + T2 * (v2 + v) / 2

                % v  = va + j*Ta*Ta/4
    % la = j/8*Ta^3 + va*Ta
% Tb = T-Ta
% lb = Tb*(v+vb)/2

%%
syms v1 v2 a j pt T1;
v  = v1 + j*T1*T1/4
T2 = 2 * sqrt((v-v2) / j )
l = T1 * (v1 + v) / 2 + T2 * (v2 + v) / 2

