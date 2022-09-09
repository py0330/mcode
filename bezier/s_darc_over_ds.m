function [darc, d2arc] = s_darc_over_ds(dp, ddp)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

darc = zeros(size(dp,2),1);
d2arc = zeros(size(dp,2),1);

for i = 1:size(dp,2)
    darc(i)  = sqrt(dp(:,i)'*dp(:,i));
    d2arc(i) = dp(:,i)'*ddp(:,i) / (darc(i));
end

end