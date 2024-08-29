%%
clear
follow_xyz = [-0.15,0.01,0.05];
target_xyz = [-0.2,0.02,0.3];
follow_v = [0,0,0];
target_v = [0.3,0.01,-0.4];
a = 12;
v = 2.5;
dt = 0.002;
%%

T = 2;

xyz_result = zeros(length(0:dt:T), 3);
target_result = zeros(length(0:dt:T), 3);
l_result = zeros(length(0:dt:T), 1);

for i = 1:length(0:dt:T)
    % 迭代更新目标位置
    last_target_xyz = target_xyz;
    target_xyz = target_xyz + target_v*dt;
%     target_xyz = target_xyz + target_v*dt + sin(i*dt/0.2)*[0.5,0.2,0.03]*dt;

    % 速度前馈
    real_target_v1 = target_v;
    if(norm(real_target_v1) > v)
        real_target_v1 = v / norm(real_target_v1) * real_target_v1;
    end

    % 补偿位置差，real_target_v2 是经过 dt 后的目标速度
    diff_xyz = target_xyz -dt*(target_v) - follow_xyz;
    l = norm(diff_xyz);
    if(l > 1e-10)
        real_target_v2 = sqrt(2*l*a) * diff_xyz / l;

        t_left = norm(real_target_v2) / a;
        if(t_left > dt)
%             target_xyz_dt = follow_xyz + real_target_v2 * dt - diff_xyz / l * a * dt *dt /2;
            real_target_v2 = real_target_v2 - diff_xyz / l * a * dt;
        else
%             target_xyz_dt = target_xyz;
            real_target_v2 = [0,0,0];
        end
    else
        real_target_v2 = [0,0,0];
    end
    if(norm(real_target_v2) > v)
        real_target_v2 = v / norm(real_target_v2) * real_target_v2;
    end
    
    % 前馈与位置补偿，构造出了目标速度
    real_target_v = real_target_v1 + real_target_v2;
    if(norm(real_target_v) > v)
        real_target_v = v / norm(real_target_v) * real_target_v;
    end
    
    % 计算加速度
    diff_v = real_target_v - follow_v;
    % 此时速度基本已经跟上目标速度，首先check 位置是否也可以满足条件
    vavg = diff_xyz/dt + real_target_v1;
    [dis, r] = s_is_in_vavg_boundage(follow_v,real_target_v,a,dt,vavg);
    if(norm(diff_v) < a * dt && dis < r)
        follow_v = real_target_v;
        follow_xyz = target_xyz;
    else
        n = norm(diff_v);

        if(n < a* dt)
            follow_a = diff_v / dt;
        else
            follow_a = diff_v / n * a;
        end

        follow_xyz = follow_xyz + follow_v*dt + follow_a * dt * dt / 2;
        follow_v = follow_v + follow_a * dt;
    end



   
%     if(n > 1e-12)
%         if(n < a* dt)
%             follow_a = v_diff / dt;
%         else
%             
%         end
%     else
%         follow_a = [0,0,0];
%     end
%     
%     if(i > 100)
%         real_target_v
%         follow_v
%         v_diff
%         111;
%     end
% 
%     next_follow_v = follow_v + follow_a * dt;
% %     next_follow_xyz = follow_xyz + (follow_v + next_follow_v)/2*dt;
%     next_follow_xyz = follow_xyz + follow_v*dt + follow_a * dt * dt / 2;
%     
%     follow_xyz = next_follow_xyz;
%     follow_v = next_follow_v;

    xyz_result(i,:) = follow_xyz;
    target_result(i,:) = target_xyz;
    l_result(i,:) = norm(diff_xyz);
end

subplot(3,1,1);
plot(xyz_result(:,:))
hold on
plot(target_result(:,:))
subplot(3,1,2);
plot(diff(xyz_result(:,:))/dt)
hold on
plot(diff(target_result(:,:))/dt)
subplot(3,1,3);
plot(diff(diff(xyz_result(:,:)))/dt/dt)
hold on
% plot(diff(diff(target_result(:,:)))/dt/dt)

vel = diff(xyz_result(:,:))/dt
acc = diff(diff(xyz_result(:,:)))/dt/dt

%%
% plot(sum(vel'.*vel'))

