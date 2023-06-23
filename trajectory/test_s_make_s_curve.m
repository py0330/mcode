%%

%% test l2
clear
a  = 3;
j  = 2;
va = 0.5;
max_v = 100;
T = 0.99;
pa=0;
pb = 0.3;
% Ta = 3;
% v  = s_acc_vend(va,a,j,Ta);
% la = Ta * (va + v) / 2;
% Tb = 3 - 1e-10;
% vb = s_acc_vend(v,-a,-j,Tb);
% lb = Tb * (v + vb) / 2;
% l  = la + lb;
% 
% pa = 1.8;
% pb = pa + l;
% T  = Ta + Tb;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,max_v,a,j,T);
%% vb不为0，达不到max_v，a段不可达 max_a，b段不可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

max_v = 100;
Ta = 3 - 1e-10;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 - 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
l  = la + lb;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,max_v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end


%% vb不为0，达不到max_v，a段不可达 max_a，b段可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

max_v = 100;
Ta = 3 - 1e-10;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 + 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
l  = la + lb;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,max_v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end


















%% vb不为0，达不到max_v，a段可达 max_a，b段不可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

max_v = 100;
Ta = 3 + 1e-10;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 - 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
l  = la + lb;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,max_v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end


















%% vb不为0，达不到max_v，a段可达 max_a，b段可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

max_v = 100;
Ta = 3 + 1e-10;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 + 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
l  = la + lb;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,max_v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end


















%% vb不为0，可达max_v，b段可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

Ta = 3.5;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 + 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
Tc = 0.5;
lc = Tc*v;
l  = la + lb + lc;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb + Tc;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end


















%% vb不为0，可达max_v，b段可达max_a
clear;
a  = 3;
j  = 2;
va = 0.5;

Ta = 3.5;
v  = s_acc_vend(va,a,j,Ta);
la = Ta * (va + v) / 2;
Tb = 3 - 1e-10;
vb = s_acc_vend(v,-a,-j,Tb);
lb = Tb * (v + vb) / 2;
Tc = 0.5;
lc = Tc*v;
l  = la + lb + lc;

pa = 1.8;
pb = pa + l;
T  = Ta + Tb + Tc;

[vb_r,v_r,Ta_r,Tb_r] = s_make_s_curve(pa,va,pb,v,a,j,T);

if(abs(vb - vb_r) > 1e-10 || abs(v_r - v) > 1e-10 ...
        || abs(Ta_r - Ta) > 1e-10|| abs(Tb_r - Tb) > 1e-10)
    error('test failed')
end

%%
pb = 0.16625342062919118
vc_max = 0.094677962311968483
vb_max = 0.080968858540273123
a = 0.25000000000000000
j = 2.5
pa = 0
va = 0
T =2.0060669712378374
vb = 0.080968858540273123
vc = 0.094677962311968455
Ta = 0.47871184924787380
Tb =0.14810322763097547
mode = 0
%%
i = 1;
dt = 0.001;
p = zeros(size(0:0.001:T))
for t=0:0.001:T
    p(i) = s_s_curve(t, pa, va, pb, vb, vc, a, j, T, Ta, Tb, mode);
    i = i+1;
end

t = 0:0.001:T
plot(t,p)

plot(t(1:end-1),diff(p)/dt)

plot(t(1:end-2),diff(diff(p))/dt/dt)


center = [0.57345962137149420;0.22504954769300603;0.23562075307425620]
p0 = [0.60198384418161333;0.20742573893289451;0.24832503880025239]
axis = [0.0019395681385318642;-0.58269493919875459;-0.81268865250327915]
radius = 0.035855666190186097

arc_at = p/radius;

rx = p0 - center;
ry = cross(axis, rx);
s = sin(arc_at);
c = cos(arc_at);
xyz = center + s .* ry + c .* rx

tool_pm = [eye(3),[0.020291056968836301;-0.0089806453668494404;0.047480744137520500]
    0,0,0,1];

maki_pm = [-0.9997534921377772   0.0000529070217595   -0.0222025259636535   0.60254987856721909
0.0000542150089666   0.9999999968303550   -0.0000583097144419   0.26249003936649468
0.0222025228082860   -0.0000594990507831   -0.9997534918372681   0.26225273101178792
0.0000000000000000   0.0000000000000000   0.0000000000000000   1.0000000000000000];

answ = (maki_pm) * [eye(3),xyz(:,1861);[0,0,0,1]] * inv(maki_pm)


diff_ = 0.020291056968836301

xyz(1,1861)
xyz(1,1861) + diff_
