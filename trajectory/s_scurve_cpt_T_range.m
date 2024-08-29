function [Tmax, Tmin] = s_scurve_cpt_T_range(pa, pb, va_upper, va_below, vb_max, vc_max, a, j)
% 计算当前点位所需的最大最小时间
%
% pa     : init pos
% va     : init vel
% pb     : end  pos
% max_vb : max  end vel
% v      : max  vel  during period
% a      : max  acc  during period
% j      : max  jerk during period
% T      : period
%
% Tmax：开始时尽可能快的减速，若减速到0，则为inf，否则以到达pb的时间为准
% Tmin：开始时尽可能快的加速，直到速度最大，之后保持最大速度到终点


% 确保必然可以实现 T
pt = pb - pa;
if(pt < (a*(a^2/j + 2*vb_max))/j)
    % clear
    % syms T va pt j T a
    % vb = va + j*T*T/4;
    % l  = T*(vb+va)/2;
    % collect(l, T)
    % solve(l==pt,T)
    T = newton_raphson_binary_search(@(T)(...
        j*T^3 + 8*vb_max*T - 8*pt) ...
        ,0,2*a/j...
        ,10*eps);
    va_upper = min(va_upper, vb_max + j*T*T/4);
else
    % clear
    % syms T va pt j T a
    % vb = va + T*a - a^2/j;
    % l  = T*(vb+va)/2;
    % collect(l, T)
    % solve(l==pt,T)
    T = (a^2 - 2*j*vb_max + 2*(a^4/4 - a^2*j*vb_max + 2*pt*a*j^2 + j^2*vb_max^2)^(1/2))/(2*a*j);
    va_upper = min(va_upper, vb_max + T*a - a^2/j);
end


Tmax = s_scurve_cpt_Tmax(pa, va_below, pb, vb_max, vc_max, a, j);
Tmin = s_scurve_cpt_Tmin(pa, va_upper, pb, vb_max, vc_max, a, j);

if(Tmin == -1)
    Tmin = s_scurve_cpt_Tmin(pa, va_upper, pb, vb_max, vc_max, a, j);
    error('error')
end

end
