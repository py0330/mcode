function [T, Ta, Tb, vb, vc, mode] = s_make_s_curve_multiple(p0, v0, pos, vb_max, vc_max, acc, jerk)
% 计算当前点位所需的最大最小时间
%
% pa        : init pos
% va        : init vel
% pb        : end  pos
% max_vend  : max end vel
% max_vel   : max vel  during period
% acc       : max acc  during period
% jerk      : max jerk during period
%
% T         : time
% Ta        : time acc
% Tb        : time dec
% real_vend : real velocity at end of each path 
% real_vel  : real max velocity of each path

[m,n] = size(pos);

Tmax      = zeros(m,n);
Tmin      = zeros(m,n);
T         = zeros(m,1);
vb        = zeros(m,n);
vb_upper  = zeros(m,n);
vb_below  = zeros(m,n);
vc        = zeros(m,n);
Ta        = zeros(m,n);
Tb        = zeros(m,n);
mode      = zeros(m,n);

% 正向迭代，计算每个节点的T
for i = 1:m
%     i
    for j = 1:n
%         j

%         if(i==5 && j == 4)
%             i
%         end

        if(i == 1)
            pa = p0(1,j);
            pas = p0;
            va_uppers = v0;
            va_belows = v0;
            va_upper = v0(1,j);
            va_below = v0(1,j);
        else
            pa = pos(i-1,j);
            pas = pos(i-1,:);
            va_uppers = vb_upper(i-1,:);
            va_belows = vb_below(i-1,:);
            va_upper = vb_upper(i-1,j);
            va_below = vb_below(i-1,j);
        end
        [Tmax(i,j), Tmin(i,j)] = s_scurve_cpt_T_range( ...
            pa, ...
            pos(i,j), ...
            va_upper, ...
            va_below, ...
            vb_max(i,j), ...
            vc_max(i,j), ...
            acc(i,j), ...
            jerk(i,j));
    end
%     pa
%     va
%     pos
%     max_vend
%     max_vel
%     acc
%     jerk
%     Tmax
%     Tmin
    

    
    Tmax_all = min(Tmax(i,:));
    Tmin_all = max(Tmin(i,:));
    
    % 二分法 search max T
    if(Tmax_all == inf)
        T_upper = 1000;
    else
        T_upper = Tmax_all;
    end
    T_below = Tmin_all;
    
    diff      = abs(T_upper - T_below);
    diff_last = diff * 2;
    while(diff < diff_last)
        diff_last = diff;

        T_next = (T_upper + T_below)/2;
        
        if(s_test_curve_slow(i, pas, va_uppers, va_belows, pos, vb_max, vc_max, acc, jerk, T_next))
            T_upper = T_next;
        else
            T_below = T_next;
        end

        diff = abs(T_upper - T_below);
    end
    T(i) = T_upper;
    

    for j=1:n
%         if(i==5 && j==5)
%             i;
%         end

        if(i > 1)
            [vb_upper(i,j), vb_below(i,j)] = s_scurve_cpt_vb_range(...
                pos(i-1,j), ...
                pos(i,j), ...
                vb_upper(i-1,j), ...
                vb_below(i-1,j), ...
                vc_max(i,j), ...
                vb_max(i,j), ...
                acc(i,j), ...
                jerk(i,j), ...
                T(i));
        else
            [vb_upper(i,j), vb_below(i,j)] = s_scurve_cpt_vb_range(...
                pas(j), ...
                pos(i,j), ...
                v0(j), ...
                v0(j), ...
                vc_max(i,j), ...
                vb_max(i,j), ...
                acc(i,j), ...
                jerk(i,j), ...
                T(i));
        end
    end
end

% 逆向迭代，计算每个节点的Ta, Tb, vc, va等信息
for i = m:-1:1
%     i
    for j = 1:n
%         j
        if(i == 1)
            pa = p0(1,j);
            va = v0(1,j);
            
            va_upper = va;
            va_below = va;

            [va, vc(i,j), Ta(i,j), Tb(i,j), mode(i,j)] = s_scurve_cpt_vavc(pa, pos(i,j), vb(i,j), va_upper, va_below, vc_max(i,j), acc(i,j), ...
                jerk(i,j), T(i));
            
%             if(vc(i,j) > max(va, vb(i,j)) || vc(i,j) < min(va, vb(i,j)))
%                 mode(i,j) = 0;
%                 Ta(i,j)=s_acc_time(va,vc(i,j),acc(i,j),jerk(i,j));
%                 Tb(i,j)=s_acc_time(vc(i,j),vb(i,j),acc(i,j),jerk(i,j));
%             else
%                 mode(i,j) = 1;
%                 if(abs(va - vb(i,j))<1e-10)
%                     Ta(i,j)=T(i)/2;
%                     Tb(i,j)=T(i)/2;
%                 else
%                     T_va_to_vb = s_acc_time(va,vb(i,j),acc(i,j),jerk(i,j));
%                     Ta(i,j)=abs((vc(i,j) - vb(i,j))/(va - vb(i,j))) * (T(i)-T_va_to_vb);
%                     Tb(i,j)=abs((vc(i,j) - va)/(va - vb(i,j))) * (T(i)-T_va_to_vb);
%                 end
%             end
        else
            pa = pos(i-1,j);
            
            if(i==8 && j == 4)
%                 i
            end

            va_upper = vb_upper(i-1,j);
            va_below = vb_below(i-1,j);
            [vb(i-1,j), vc(i,j), Ta(i,j), Tb(i,j), mode(i,j)] = s_scurve_cpt_vavc(pa, pos(i,j), vb(i,j), va_upper, va_below, vc_max(i,j), acc(i,j), ...
                jerk(i,j), T(i));

%             if(vc(i,j) > max(vb(i-1,j), vb(i,j)) || vc(i,j) < min(vb(i-1,j), vb(i,j)))
%                 mode(i,j) = 0;
%                 Ta(i,j)=s_acc_time(vb(i-1,j),vc(i,j),acc(i,j),jerk(i,j));
%                 Tb(i,j)=s_acc_time(vb(i,j),vc(i,j),acc(i,j),jerk(i,j));
%             else
%                 mode(i,j) = 1;
%                 if(abs(vb(i-1,j) - vb(i,j))<10*eps)
%                     Ta(i,j)=T(i)/2;
%                     Tb(i,j)=T(i)/2;
%                 else
%                     T_va_to_vb = s_acc_time(vb(i-1,j),vb(i,j),acc(i,j),jerk(i,j));
%                     Ta(i,j)=abs((vc(i,j) - vb(i,j))/(vb(i-1,j) - vb(i,j))) * (T(i)-T_va_to_vb);
%                     Tb(i,j)=abs((vc(i,j) - vb(i-1,j))/(vb(i-1,j) - vb(i,j))) * (T(i)-T_va_to_vb);
%                 end
%             end
        end


    end
end


end

