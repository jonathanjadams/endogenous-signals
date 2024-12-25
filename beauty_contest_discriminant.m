function D = beauty_contest_discriminant(tau_u,tau_v,alpha,varphi)

% This function calculates the cubic discriminant in the beauty contest indeterminacy model as a component of MMIIES

% Date: 11/17/2024
% Contact: adamsjonathanj@gmail.com

% general cubic: a x^3 + b x^2 + c x + d = 0
a = tau_u*(1-alpha);
b = -tau_u*varphi;
c = (1+tau_v*(1-alpha));
d = -varphi*(1+tau_v);

D = b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d;

%D = varphi^4*tau_u^2*(1+tau_v)^2 + 18*tau_u^2*(1+tau_v)*(1+(1-alpha)*tau_v)*(1-alpha)*varphi^2 - 4*varphi^4*tau_u^3*(1+tau_v) - 4*tau_u*(1+(1-alpha)*tau_v)^3*(1-alpha)-27*tau_u^2*(1+tau_v)^2*(1-alpha)^2*varphi^2;
%why did I code it like that?  The discriminant was wrong

%Wolfram GPT confirms my new solution
end

