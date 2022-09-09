function t = s_acc_time(va,vb,a,j)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
    % CASE 1 |v - va| >  a^2 / j, 此时加速段可以达到最大加速度：
    %        T0 = |v - va|/a + a/j            |    |v - va| >  a^2 / j
    % CASE 2 |v - va| <= a^2 / j, 此时加速段（或减速段）无法达到最大加速度：
    %        T1 = 2 * sqrt( |v - va| / j )    |    |v - va| <= a^2 / j

    if(abs(vb-va)>a^2 / j)
        t = abs(vb-va) / a + a/j;
    else
        t = 2 * sqrt( abs(vb-va) / j );
    end
end
