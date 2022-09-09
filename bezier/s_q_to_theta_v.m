function [theta,v,dtheta,dv,ddtheta,ddv] = s_q_to_theta_v(q,dq,ddq)
%UNTITLED6 此处提供此函数的摘要
%   此处提供详细说明

q=q(:);
dq=dq(:);
if(q(4) < 0)
    q = -q;
    dq = -dq;
end

theta = acos(q(4));
if(theta > 1e-7)
    v = q(1:3)/sin(theta);
    
    dtheta = -dq(4)/sin(theta);
    dv = (dq(1:3)-v*cos(theta)*dtheta)/sin(theta);
    
    ddtheta = -(ddq(4)*sin(theta) - dq(4)*cos(theta)*dtheta)/sin(theta)^2;
    ddv = (...
            (ddq(1:3)-dv*cos(theta)*dtheta + v*sin(theta)*dtheta^2-v*cos(theta)*ddtheta)*sin(theta)...
           -(dq(1:3)-v*cos(theta)*dtheta)*cos(theta)*dtheta...
          )/sin(theta)^2;
else
    % 这里不对，后面还需要再推导
    v = dq(1:3)/norm(dq(1:3));
    
    dtheta = (v'*dq(1:3))/cos(theta);
    dv = zeros(3,1);
    
    ddtheta = 0;
    ddv = zeros(3,1);

end





end