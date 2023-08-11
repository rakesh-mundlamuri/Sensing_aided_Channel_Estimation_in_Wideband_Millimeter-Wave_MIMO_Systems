function [rho_rc] = raised_cosine_filter(d,delay,Ts,beta)

% if(abs(d*Ts-delay) == Ts/(2*beta))
%     rho_rc = (pi/(4*Ts))*sinc(1/(2*beta));
%     %rho_rc = (pi/(4))*sin(1/(2*beta))/(1/(2*beta));
% else
    rho_rc = sinc((d*Ts-delay)/Ts).*(cos(pi*beta*(d*Ts-delay)/Ts)./(1-((2*beta*(d*Ts-delay)/Ts).^2)));
    %rho_rc = (sin((d*Ts-delay)/Ts)/((d*Ts-delay)/Ts))*(cos(pi*beta*(d*Ts-delay)/Ts)/(1-((2*beta*(d*Ts-delay)/Ts)^2)));
 % end

end
