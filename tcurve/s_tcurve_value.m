function p_at = s_tcurve_value(pb, pe, vb, ve, T, Ta, Tb, mode, v, a, t)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

if(mode == 0)
    if(t < Ta)
        p_at = pb + vb*t + a*t.*t/2;
    elseif(t < T - Tb)
        p_at = pb + vb*Ta + a*Ta*Ta/2 + v*(t-Ta);
    else
        p_at = pe - ve*(T-t) - a*(T-t).*(T-t)/2;
    end
else
    if(t < Ta)
        p_at = pb + vb*t;
    elseif(t < T - Tb)
        p_at = pb + vb*t + a*(t-Ta).*(t-Ta)/2;
    else
        p_at = pe - ve*(T-t);
    end
end




end