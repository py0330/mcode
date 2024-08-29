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

Z1 = a^2/j;
Z2 = T^2*j;

T_va_vcmax = s_acc_time(va,vc_max,a,j);
T_vcmax_vbmax = s_acc_time(vc_max,vb_max,a,j);
T_vcmax_0     = s_acc_time(vc_max,0,a,j);
T_va_0     = s_acc_time(va,0,a,j);
T_0_vbmax = s_acc_time(0,vb_max,a,j);
T_acc     = 2*a/j;

l0 = inf;

%-----------------l1------------------------%
if(T>=T_va_vcmax+T_vcmax_vbmax)
    Ta = T_va_vcmax;
    Tb = T_vcmax_vbmax;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = vb_max;
    l1 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;

    if(l1 <= pt)
        vb_below = vb_max;
        return;
    end
else
    l1 = l0;
end
%-----------------l2------------------------%
if(T>=T_va_vcmax+T_acc)
    Ta = T_va_vcmax;
    Tb = 2*a/j;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = max(vc_max-a*a/j,0);

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l2 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
        if(l2 <= pt)
            la = Ta*(va+vc_max)/2;
            Tb = max((la + vc_max*(T - Ta) - pt)*8/j,0)^(1/3);
            vb = vc_max - j*Tb*Tb/4;

            vb_below = vb;
            % debug check %
            if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
                error('s_scurve_cpt_vb_below failed: l2')
            end
            return;
        end
    else
        l2 = l1;
    end
else
    l2 = l1;
end
%-----------------l3------------------------%
if(T>=T_va_vcmax+T_vcmax_0)
    Ta = T_va_vcmax;
    Tb = T_0_vcmax;
    Tc = T-Ta-Tb;
    vc = vc_max;
    vb = 0;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l3 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;

        if(l3 <= pt)
            la = Ta*(va+vc_max)/2;
            B = -a/j;
            C = -(2*la - 2*pt + 2*vc_max*(T - Ta))/a;

            Tb = max((-B+sqrt(max(B^2-4*C,0)))/2,0);
            vb = vc_max - Tb*a + Z1;

            vb_below = vb;
            % debug check %
            if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
                error('s_scurve_cpt_vb_below failed: l3')
            end
            return;
        end
    else
        l3 = l2;
    end
else
    l3 = l2;
end
%-----------------l4------------------------%
if(T>=T_va_vcmax + T_vcmax_vbmax)
    Ta = T_va_vcmax;
    Tb = T-Ta;
    Tc = 0;
    vc = vc_max;
    vb = s_acc_vend(vc_max,-a,-j,Tb);

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l4 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
        if(l4 <= pt)
            la = Ta*(va+vc_max)/2;
            B = -a/j;
            C = -(2*la - 2*pt + 2*vc_max*(T - Ta))/a;

            Tb = max((-B+sqrt(max(B^2-4*C,0)))/2,0);
            vb = vc_max - Tb*a + Z1;

            vb_below = vb;
            % debug check %
            if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
                error('s_scurve_cpt_vb_below failed: l4')
            end
            return;
        end
    else
        l4 = l3;
    end
else
    l4 = l3;
end
%-----------------l5------------------------%
if(T<=T_va_vcmax)
    Ta = T;
    Tb = 0;
    Tc = 0;
    vc = s_acc_vend(va,a,j,Ta);
    vb = vc;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l5 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
        if(l5 <= pt)
            vb = s_acc_vend(va,a,j,T);

            vb_below = vb;
            % debug check %
            if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
                error('s_scurve_cpt_vb_below failed: l5')
            end
            return;
        end

    else
        l5 = l4;
    end
else
    l5 = l4;
end
%-----------------l6------------------------%
if(T>=4*a/j)
    Ta = T-2*a/j;
    Tb = 2*a/j;
    Tc = 0;
    vc = s_acc_vend(va,a,j,Ta);
    vb = vc-a*a/j;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l6 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;

        if(l6 <= pt)
            vb = s_acc_vend(va,a,j,T);

            vb_below = vb;
            % debug check %
            if(abs(Tb*vb_below/2 + vc_max*(T-Ta-Tb) + Ta*va/2 -pt) > cons)
                error('s_scurve_cpt_vb_below failed: l5')
            end
            return;
        end
    else
        l6 = l5;
    end
else
    l6 = l5;
end
%-----------------l7------------------------%
if(T>=4*a/j)
    Ta = 2*a/j;
    Tb = T-2*a/j;
    Tc = 0;
    vc = s_acc_vend(va,a,j,Ta);
    vb = s_acc_vend(va,-a,-j,Tb);

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l7 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l7 = l6;
    end
else
    l7 = l6;
end
%-----------------l8------------------------%
if(T<=4*a/j)
    Ta = 2*a/j;
    Tb = T-2*a/j;
    Tc = 0;
    vc = s_acc_vend(va,a,j,Ta);
    vb = vc;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l8 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l8 = l7;
    end
else
    l8 = l7;
end
%-----------------l9------------------------%
if(T<=4*a/j)
    Ta = T-2*a/j;
    Tb = 2*a/j;
    Tc = 0;
    vc = s_acc_vend(va,a,j,Ta);
    vb = vc;

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l9 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l9 = l8;
    end
else
    l9 = l8;
end
%-----------------l10------------------------%
if(1)
    Ta = 0;
    Tb = T;
    Tc = 0;
    vc = va;
    vb = s_acc_vend(va,a,j,Tb);

    if(0 <= vb && vb <= vb_max && 0 <= vc && vc <= vc_max)
        l10 = (va+vc)*Ta/2 + (vb+vc)*Tb/2+vc*Tc;
    else
        l10 = l9;
    end
else
    l10 = l9;
end

end