function [theta] = thetacalcy(y,w1, Lc1, L1)
%thetacalcy Calculate theta when piston 1 position is known
%   The function uses the bracketing root finding method to calculate theta
%   from y. can also calculate from yp2, but sign of w1 should be negative.
%% Calculate theta using root finding method (bracketing)
% thetamax = 45;
% thetamin = -45;
% thetaguess = 0;
% rootval = 100;
% cutoff = 0.000005;
% y;
% while abs(rootval) > cutoff
%     rootval = -w1*sind(thetaguess) + Lc1*(1-sqrt(1-((L1-w1*cosd(thetaguess))/H1)^2)) - y;
%     if rootval < -cutoff
%         thetamax = thetaguess;
%         thetaguess = (thetamax + thetamin) /2;
%     elseif rootval > cutoff
%         thetamin = thetaguess;
%         thetaguess = (thetamax + thetamin) /2;
%     end
% end
% 
% 
% 
% theta = thetaguess;
% end

%% Calculate theta using Newton-Raphson method
cutoff = 0.0000005;
iterationcap = 2000;
theta = 0;
oldtheta = theta;
iteration = 0;
residual = cutoff+1;
while residual > cutoff && iteration < iterationcap
    
    ftheta = -w1*sind(theta) + Lc1*(1-sqrt(1-((L1-abs(w1)*cosd(theta))/Lc1)^2))-y;
    fprimetheta = (abs(w1)*sind(theta)*(L1-abs(w1)*cosd(theta)))/(Lc1*sqrt(1-((L1-abs(w1)*cosd(theta))/Lc1)^2))-w1*cosd(theta);
    oldtheta = theta;
    theta = theta - ftheta/fprimetheta;
    residual = abs((theta-oldtheta)/oldtheta);
    iteration = iteration + 1;
end