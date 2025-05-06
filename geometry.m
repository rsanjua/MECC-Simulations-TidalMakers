% geometry calls
function [y, ydot, ydotdot, theta, thetadot, thetadotdot, yp2, yp2dot, yp2dotdot,phi1,phi2] = geometry(y,ydot,yn1,yn2,yn3,yp2,yp2dot,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,known)
% this function calculates the positional variables and their derivatives
% of the entire system based on previous information and one piece of
% current information noted by the variable known
    switch known
        case "y"
            % y = y;
            ydot = (3*y - 4*yn1 + yn2)/(2*c.delt); % ydot
            theta = thetacalcy(y,c.w1, c.Lc1, c.L1); % 17a
            yp2 = c.w2*sind(theta) + c.Lc2*(1-sqrt(1-((c.L2-c.w2*cosd(theta))/c.H2)^2)); % 17b
            yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt); % yp2dot
    
        case "ydot"
            y = (2*c.delt*ydot + 4*yn1 - yn2)/3; %ydot
            % ydot = ydot;
            theta = thetacalcy(y,c.w1, c.Lc1, c.L1); % 17a
            yp2 = c.w2*sind(theta) + c.Lc2*(1-sqrt(1-((c.L2-c.w2*cosd(theta))/c.H2)^2)); % 17b
            yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt); % yp2dot
        
        case "yp2dot"
            
            yp2 = (2*c.delt*yp2dot + 4*yp2n1 - yp2n2)/3;
            theta = thetacalcy(-yp2,c.w2, c.Lc2, c.L2); % 17a
            y = -1*c.w1*sind(theta) + c.Lc1*(1-sqrt(1-((c.Lc1-c.w1*cosd(theta))/c.H1)^2));
            % yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt) % yp2dot
            ydot = (3*y - 4*yn1 + yn2)/(2*c.delt);
            
    end
   
    ydotdot = (2*y - 5*yn1 + 4*yn2 - yn3) / (c.delt^2); % ydotdot
    thetadot = (3*theta - 4*thetan1 + thetan2)/(2*c.delt); % thetadot
    yp2dotdot = (2*yp2 - 5*yp2n1 + 4*yp2n2 - yp2n3) / (c.delt^2); % yp2dotdot
    thetadotdot = (2*theta - 5*thetan1 + 4*thetan2 - thetan3) / (c.delt^2); % thetadotdot
    phi2 = asind((c.L2-c.w2*cosd(theta))/c.H2);
    phi1 = asind((c.L1-c.w1*cosd(theta))/c.H1);


end