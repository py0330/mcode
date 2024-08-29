function [Ta,Tb,v,a,mode] = s_tcurve_param(pb, pe, vb, ve, vmax, amax, T)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

pt = pe - pb;

Tb2e = abs(vb-ve)/amax;
pb2e = (vb+ve)/2*Tb2e;

if((pt-pb2e) > (T-Tb2e)*max(vb,ve))
    mode = 0;
    a   = amax;

    A = 1;
    B = -T*a - vb - ve;
    C = (vb*vb+ve*ve)/2+pt*a;

    v   = (-B-sqrt(B*B-4*A*C))/(2*A);

    Ta  = (v-vb)/a;
    Tb  = (v-ve)/a;
elseif((pt-pb2e) > (T-Tb2e)*min(vb,ve))
    mode = 1;
    v  = (pt-pb2e) / (T-Tb2e);
%     Ta = (v-ve)/(vb-ve)*(T-Tb2e);
%     Tb = (vb-v)/(vb-ve)*(T-Tb2e);
    Ta = (pt-pb2e)/(vb-ve) - ve/(vb-ve)*(T-Tb2e);
    Ta = max(0,Ta);
    Ta = min(Ta, T-Tb2e);
    Tb = T-Tb2e - Ta;
    a  = sign(ve-vb)*amax;
else
    mode = 0;
    a   = -amax;
    
    %
    % A*v^2 + B*v + C == 0
    % v^2 - (T*a+vb+ve)*v + (vb*vb+ve*ve)/2+pt*a
    % Ta = (v-vb)/a
    % Tb = (v-ve)/a
    % Ta*(v+vb)/2 + Tb*(v+ve)/2 + (T-Ta-Tb)*v - pt
    % v^2 - a*(T + vb/a + ve/a)*v + a*(pt + vb^2/(2*a) + ve^2/(2*a))
    % (-1/a)*v^2 + (T + vb/a + ve/a)*v - vb^2/(2*a) - ve^2/(2*a) - pt

    A = 1;
    B = -T*a - vb - ve;
    C = (vb*vb+ve*ve)/2+pt*a;

    v = (-B+sqrt(B*B-4*A*C))/(2*A);
    
%     v = (- vb^2/(2*a) + ve^2/(2*a) - pt)/((T + vb/a - ve/a));

    Ta  = (v-vb)/a;
    Tb  = (v-ve)/a;
end



end