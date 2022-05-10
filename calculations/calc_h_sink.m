function h = calc_h_sink(rho, k, mu, Pr, beta, D_h, Th_in, Th_out, Tc)
% Calculate heat transfer coefficient of heatsink fluid (eg. sea) for
% horizontal cylindrical pipe, via natural convection
g    = 9.81;				% Gravitational acceleration [m/s^2]
Ts_h = (Th_in + Th_out)/2;	% Avg. pipe surface temperature [K]

% Rayleigh number
Ra = ( (g*beta*(Ts_h - Tc)*D_h^3)/((mu/rho)^2) )*Pr;

% Nusselt correl. for nat. convect. of horizontal cylinder [Churchill and Chu]
Nu = (  0.6  +  ( 0.387*Ra^(1/6) ) / ( (1 + (0.559/Pr)^(9/16))^(8/27) )  )^2;

% Heat transfer coefficient [W/(m^2 K)]
h = k*Nu/D_h;  
end