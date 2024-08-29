function [T1, T2, T3] = s_tcurve_T_range(pb, pe, vb, ve, vmax, amax)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

% 可行的T 区间位于(T1,T2) 并 (T3,inf) 
% T1 是最短可能的时间
% T3 是之后必定可以加速到0的时间，因此可以在0速度处停留无穷的时间
% T2 是不加速到0，最长可能的时间
% 
% 【1】求 T1
%
% pt = Ta*(vb+v)/2 + Tb*(ve+v)/2 + (T-Ta-Tb)*v
%    = Ta*(vb-v)/2 + Tb*(ve-v)/2 + T*v
% 
% Ta = abs(v-vb)/amax
% Tb = abs(v-ve)/amax
%
% COND1 若(vb+ve)/2 * abs(vb-ve)/amax < pt
%    此时 v > max(vb,ve)
%    SUB 1 (vmax+vb)/2*(vmax-vb)/amax + (vmax+ve)/2*(vmax-ve)/amax < pt
%       v  = vmax
%       T1 = (pt - Ta*(vb+v)/2 - Tb*(ve+v)/2)/v
%    SUB 2 否则
%       v^2/amax - vb^2/(2*amax) - ve^2/(2*amax) == pt
%       sig = sign(pt-va/2*abs(va)/amax-vb/2*abs(vb)/amax)
%       v  = sig*sqrt(pt*amax + vb^2/2 + ve^2/2)
%       T1 = (v-vb)/amax + (v-ve)/amax
% COND2 否则
%    v < min(vb+ve)
%    SUB 1 (-vmax+vb)/2*(-vmax-vb)/amax + (-vmax+ve)/2*(-vmax-ve)/amax > pt
%       v  = -vmax
%       T1 = (pt - Ta*(vb+v)/2 - Tb*(ve+v)/2)/v
%    SUB 2 否则
%       v^2/amax - vb^2/(2*amax) - ve^2/(2*amax) == pt
%       sig = sign(pt-va/2*abs(va)/amax-vb/2*abs(vb)/amax)
%       v  = sig*sqrt(-pt*amax + vb^2/2 + ve^2/2)
%       T1 = (vb-v)/amax + (ve-v)/amax
%
% 【2】求 T2
% 在 va vb v同号，且 abs(v) > max(abs(vb), abs(ve))，可能会出现 T2 
% COND1 va vb v > 0 AND v > max(vb,ve)
%    sig = pt+va/2*abs(va)/amax+vb/2*abs(vb)/amax
%    若 sig > 0
%    v  = sig*sqrt(pt*amax + vb^2/2 + ve^2/2)
%    T2 = 
%
%
%

pt = pe - pb;

Tb2e = abs(vb-ve)/amax;
pb2e = (vb+ve)/2*Tb2e;

if(pb2e < pt)
    if((vmax+vb)/2*(vmax-vb)/amax + (vmax+ve)/2*(vmax-ve)/amax < pt)
        v  = vmax;
        Ta = (v-vb)/amax;
        Tb = (v-ve)/amax;
        T1 = (pt - Ta*(vb+v)/2 - Tb*(ve+v)/2)/v + Ta + Tb;
    else
        v  = sqrt(pt*amax + vb^2/2 + ve^2/2);
        T1 = (v-vb)/amax + (v-ve)/amax;
    end
else
    if((-vmax+vb)/2*(vb+vmax)/amax + (-vmax+ve)/2*(ve+vmax)/amax > pt)
        v  = -vmax;
        Ta = (vb-v)/amax;
        Tb = (ve-v)/amax;
        T1 = (pt - Ta*(vb+v)/2 - Tb*(ve+v)/2)/v + Ta + Tb;
    else
        v  = -sqrt(-pt*amax + vb^2/2 + ve^2/2);
        T1 = (vb-v)/amax + (ve-v)/amax;
    end
end

T2 = T1;
T3 = T1;

if(v>0 && vb>0 && ve>0)
    if(vb*vb/2/amax + ve*ve/2/amax > pt)
        v  = sqrt(-pt*amax + vb^2/2 + ve^2/2);
        T2 = (vb-v)/amax + (ve-v)/amax;
        T3 = (vb+v)/amax + (ve+v)/amax;
    end
elseif(v<0 && vb<0 && ve<0)
    if(-vb*vb/2/amax - ve*ve/2/amax < pt)
        v  = -sqrt(pt*amax + vb^2/2 + ve^2/2);
        T2 = abs(vb-v)/amax + abs(ve-v)/amax;
        T3 = abs(vb+v)/amax + abs(ve+v)/amax;
    end
end





end