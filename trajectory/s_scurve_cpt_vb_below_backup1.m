function [vb_below] = s_scurve_cpt_vb_below(pa, va, pb, vc_max, vb_max, a, j, T)
% 在给定过程中起始速度va, 最大速度v，加速度a，跃度j，时间长度T的情况下
% 自动计算末端所可能达到的【最大的】vb
%
% t      : current time
% pa     : init pos
% va     : init vel
% pb     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period
%
cons = 100*eps;
pt = pb - pa;

T_va_vcmax = s_acc_time(va,vc_max,a,j);
T_vcmax_vbmax = s_acc_time(vc_max,vb_max,a,j);
T_vcmax_0     = s_acc_time(vc_max,0,a,j);
T_va_0     = s_acc_time(va,0,a,j);
T_0_vbmax = s_acc_time(0,vb_max,a,j);
T_acc     = 2*a/j;

l0 = inf;

if(T>=T_va_vcmax+T_vcmax_vbmax)
    Ta = T_va_vcmax;
    Tb = T_vcmax_vbmax;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = vb_max;
    l1 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
else
    l1 = l0;
end

if(T>=T_va_vcmax+T_acc)
    Ta = T_va_vcmax;
    Tb = 2*a/j;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = vc_max-a*a/j;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l2 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l2 = l1;
    end
else
    l2 = l1;
end

if(T>=T_va_vcmax+T_0_vcmax)
    Ta = T_va_vcmax;
    Tb = T_0_vcmax;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = 0;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l3 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l3 = l2;
    end
else
    l3 = l2;
end

if(T>=T_va_vcmax)
    Ta = T_va_vcmax;
    Tb = T-Ta;
    Tc = 0;
    vc = vc_max;
    vb = s_acc_vend(vc_max,-a,-j,Tb);

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l4 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l4 = l3;
    end
else
    l4 = l3;
end

end