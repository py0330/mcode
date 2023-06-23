function [ds, d2s, darc, d2arc] = s_ds_over_darc(dp, ddp)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

ds = zeros(size(dp,2),1);
d2s = zeros(size(dp,2),1);
darc = zeros(size(dp,2),1);
d2arc = zeros(size(dp,2),1);

for i = 1:size(dp,2)
    % 
    % ds / dA = 1 / (dA / ds)
    % d2s / dA2 = - [d(dA/ds) / dA]  /  (dA / ds)^2
    %           = - [d2A/ds2 * ds/dA] / (dA / ds)^2
    %           = - (d2A/ds2) / (dA/ds)^3

    darc(i)  = sqrt(dp(1,i)^2 + dp(2,i)^2 + dp(3,i)^2);
    ds(i)  = 1/darc(i);

    if(darc(i) < eps)
        d2arc(i) = 0;
        d2s(i) = 0;
    else
        d2arc(i) = dp(:,i)'*ddp(:,i) / (darc(i));
        d2s(i) = -d2arc(i)./darc(i).^3;
    end
    
end

end