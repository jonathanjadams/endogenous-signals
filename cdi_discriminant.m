function D = cdi_discriminant(tau_y,tau_v,beta)

% This function calculates the cubic discriminant in the idiosycnratic confounding dynamics model as a component of MMIIES

% Date: 11/17/2024
% Contact: adamsjonathanj@gmail.com

% general cubic: a x^3 + b x^2 + c x + d = 0
a = tau_v;
b = 0;
c = tau_y+1-beta*tau_v;
d = -beta*tau_y;

%D = -4 \tau_v(\tau_y + 1 - \beta \tau_v)^3 - 27\tau_v^2 \beta^2 \tau_y^2 
D = b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d;

end

