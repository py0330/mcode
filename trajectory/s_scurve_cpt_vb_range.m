function [vb_upper, vb_below] = s_scurve_cpt_vb_range(pa, pb, va_upper, va_below, vc_max, vb_max, a, j, T)

pt = pb - pa;
vc_min = 0;
vb_min = 0;

va_upper_ori = va_upper;
va_below_ori = va_below;

% 计算恰好可以到达 pt T 的 v1 v2
v_diff = s_acc_vend(0,a,j,T);
v_dis  = pt/T - v_diff/2;

v1_below = v_dis;
v2_upper = v_diff + v_dis;

% 整个过程的速度应该在上述 v1 v2 之间，因此调整 vc 的范围
vc_min = max(vc_min, v1_below);
vc_max = min(vc_max, v2_upper);

% 但是上述 vc_min, vc_max 仍然未必可达，还需再考虑限制
v_avg = pt/T;
if(abs(vc_min - v_avg) < abs(vc_max - v_avg))
    % vc_max 可能无法达到，vc_min 一定可以达到
    vc_max = s_scurve_cpt_vb_upper(pa, vc_min, pb, vc_max, vc_max, a, j, T);
else
    % vc_min 可能无法达到
    vc_min = s_scurve_cpt_vb_below(pa, vc_max, pb, vc_max, vc_max, a, j, T);
end

% 修正 va_range
va_upper = min(va_upper, vc_max);
va_below = max(va_below, vc_min);

% 计算 vb_range
vb_upper = s_scurve_cpt_vb_upper(pa, va_below, pb, vc_max, vb_max, a, j, T);
vb_below = s_scurve_cpt_vb_below(pa, va_upper, pb, vc_max, vb_max, a, j, T);

if(vb_below > vb_upper + 1e-10)
    vb_upper = s_scurve_cpt_vb_upper(pa, va_below_ori, pb, vc_max, vb_max, a, j, T);
    vb_below = s_scurve_cpt_vb_below(pa, va_upper_ori, pb, vc_max, vb_max, a, j, T);
    error('failed in vb range');
end

% 根据 vb_max 进行修正
vb_upper = min(vb_max, vb_upper);
vb_upper = max(vb_min, vb_upper);
vb_below = max(vb_min, vb_below);
vb_below = min(vb_max, vb_below);

return;

end