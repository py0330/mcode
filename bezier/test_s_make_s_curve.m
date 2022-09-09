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


