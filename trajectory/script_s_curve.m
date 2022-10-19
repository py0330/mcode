%%
T = 3.87953494659147;
% T = 3.0;
v = 5.1;
a = 1;
j = 1;
va = 2.1;
pa = 0.8;
pb = 45;
max_vb = 3;

[Tmax, Tmin] = s_s_curve_Tmax_Tmin(pa, va, pb, max_vb, v, a, j)

if(Tmax~=inf)
    T = Tmax;
else
    T = 50;
end

T = Tmin;

dt = 0.01;
t = 0:dt:T;

p= zeros(size(t));
v1 = 0;
[vb,real_v,Ta,Tb]=s_make_s_curve(pa,va,pb,v,a,j,T);

for i=1:length(t)
[p(i),v1] = s_s_curve(t(i), pa, va, Ta, pb, vb, Tb, real_v, a, j, T);
end

Ta
Tb

v = diff(p)/dt;
a = diff(v)/dt;

subplot(1,3,1);
plot(t,p);
axis equal

subplot(1,3,2);
plot(t(2:end)-dt/2,v);
axis equal

subplot(1,3,3);
plot(t(3:end)-dt/2*3,a);
axis equal




