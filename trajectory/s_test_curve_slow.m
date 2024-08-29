function success = s_test_curve_slow(i, p0, va_upper, va_below, pos, vb_max, vel, acc, jerk, T)
% 尝试将速度降下来
%
% pa     : current pos
% va     : current vel
% pos    : target positions
% max_vb : max vel at target positions
% vel    : max vel  during period
% acc    : max acc  during period
% jerk   : max jerk during period

[m,n] = size(pos);

Tmax   = zeros(1,n);
Tmin   = zeros(1,n);

if(T == inf)
    T=1000;
end

% 因为输入的 T 对本次规划肯定是成功的，因此如果是最后一个周期，那么就直接成功
if(i == m)
    success = 1;
    return;
else
    for j=1:n
%         'in test'
%         i
        % 计算 vb_upper vb_below，他们也是下一次的 va，因此放到同一组变量中
        va_upper_ori = va_upper;
        va_below_ori = va_below;
        [va_upper(j), va_below(j)] = s_scurve_cpt_vb_range(...
            p0(j), ...
            pos(i,j), ...
            va_upper(j), ...
            va_below(j), ...
            vel(i,j),...
            vb_max(i,j),...
            acc(i,j),...
            jerk(i,j),...
            T);
        
        if(va_below(j) > va_upper(j) + 1e-10)
            [va_upper(j), va_below(j)] = s_scurve_cpt_vb_range(...
            p0(j), ...
            pos(i,j), ...
            va_upper_ori(j), ...
            va_below_ori(j), ...
            vel(i,j),...
            vb_max(i,j),...
            acc(i,j),...
            jerk(i,j),...
            T);


            error('failed: T should make vb range success');
            success = 0;
            return;
        end

        p0(j) = pos(i,j);

        [Tmax(j), Tmin(j)] = s_scurve_cpt_T_range( ...
            p0(j), ...
            pos(i+1,j), ...
            va_upper(j), ...
            va_below(j), ...
            vb_max(i+1,j), ...
            vel(i+1,j), ...
            acc(i+1,j), ...
            jerk(i+1,j));
    end

    Tmax_all = min(Tmax(:));
    Tmin_all = max(Tmin(:));

    if(Tmax_all == inf)
%         fprintf('test slow success:%d \n T:%f \n Tmax_ALL:%f \n Tmin_all:%f',i,T,Tmax_all,Tmin_all)
%         Tmax
%         Tmin
%         v0
        success = 1;
        return;
    end

    if(Tmax_all == -1 || Tmax_all < Tmin_all)
%         fprintf('test slow failed:%d \n T:%f \n Tmax_ALL:%f \n Tmin_all:%f',i,T,Tmax_all,Tmin_all)
%         Tmax
%         Tmin
%         v0
        success = 0;
        return;
    end
    
    success = s_test_curve_slow(i + 1, p0, va_upper, va_below, pos, vb_max, vel, acc, jerk, Tmax_all);
    return;

end



end


