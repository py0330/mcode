function [p, v, a, j] = s_s_curve(t, pa, va, pb, vb, v, a, j, T, Ta, Tb, mode)
% 在给定过程中最大速度v，加速度a，跃度j，时间长度T的情况下
% 自动计算在最小的末端速度v1的情况下，整个过程的s曲线
%
% t      : current time
% pa     : init pos
% va     : init vel
% pb     : end  pos
% v      : max vel  during period
% a      : max acc  during period
% j      : max jerk during period
% T      : period

%CASE B
if(mode == 1)
    if(t < Ta)
        p = pa + va * t;
    elseif(t < T - Tb)
        Tacc = T - Ta - Tb;
        si = sign(vb - va);
        t = t - Ta;
        if(Tacc >= 2*a/j)
            if(t < a/j)
                p = pa + va * Ta + va * t + si/6*j*t^3;
            elseif(t < (Tacc - a/j))
                p = pa + va * Ta + va * a / j + si/6 * a^3 / j^2 + (va + si/2*a^2/j) * (t - a/j)+ si*a/2*(t-a/j)^2;
            else
                p = pa + va * Ta + (va+vb)/2*Tacc-(vb*(Tacc-t) - si/6*j*(Tacc-t)^3);
            end
        else
            if(t < Tacc / 2)
                p = pa + va * Ta + va * t + si/6*j*t^3;
            else
                p = pa + va * Ta + (va+vb)/2*Tacc - (vb * (Tacc-t) - si/6*j*(Tacc-t)^3);
            end
        end
    else
        p = pb - vb * (T-t);
    end
    return;
end

%CASE A
if(t < Ta)
    si = sign(v - va);

    if(Ta >= 2*a/j)
        if(t < a/j)
            p = pa + va * t + si/6*j*t^3;
        elseif(t < (Ta - a/j))
            p = pa + va * a / j + si/6 * a^3 / j^2 + (va + si/2*a^2/j) * (t - a/j)+ si*a/2*(t-a/j)^2; 
        else
            p = pa + (va+v)/2*Ta-(v*(Ta-t) - si/6*j*(Ta-t)^3);
        end
    else
        if(t < Ta / 2)
            p = pa + va * t + si/6*j*t^3;
        else
            p = pa + (va+v)/2*Ta - (v * (Ta-t) - si/6*j*(Ta-t)^3);
        end
    end
elseif(t < T-Tb)
    p = pa + (va+v)/2*Ta + v * (t - Ta);
else
    si = sign(vb - v);
    if(Tb >= 2*a/j)
        if(T-t < a/j)
            p = pb - vb * (T-t) + si/6*j*(T-t)^3;
        elseif(T-t < (Tb - a/j))
            p = pb - vb * a / j + si/6 * a^3 / j^2 - (vb - si/2*a^2/j) * (T - t - a/j) + si*a/2*(T - t - a/j)^2;
        else
            p = pb - (v+vb)/2*Tb + (v*(Tb-T+t) + si/6*j*(Tb-T+t)^3);
        end
    else
        if(T-t < Tb / 2)
            p = pb - vb * (T-t) + si/6*j*(T-t)^3;
        else
            p = pb - (v+vb)/2*Tb + (v*(Tb-T+t) + si/6*j*(Tb-T+t)^3);
        end
    end
end

end

