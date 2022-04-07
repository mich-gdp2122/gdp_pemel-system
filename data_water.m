function [rho, cp, k, mu, Pr] = data_water(T_in, T_out)
% Properties of water, averaged b/tw 0.01 & 0.1 MPa

% Rows of the tables correspond to the temperature vector
%%
T = mean([T_in T_out]);	% Convert to avg temperature
T_TLU = 273.1600:10:373.16; % Temperature vector [K]

rho_TLU = [ ...
    mean([999.8    999.8]) ...
    mean([999.7    999.7]) ...
    mean([998.2    998.2]) ...
    mean([995.6    995.6]) ...
    mean([992.2    992.2]) ...
    mean([988.0    988.0]) ...
    mean([983.2    983.2]) ...
    mean([977.8    977.8]) ...
    mean([971.8    971.8]) ...
    mean([965.3    965.3]) ...
    mean([958.8    958.8]) ...
    ]; % Density table [kg/m^3]

cp_TLU  = 1000*[ ...
    mean([4.2199    4.2194]) ...
    mean([4.1955    4.1951]) ...
    mean([4.1843    4.1840]) ...
    mean([4.1801    4.1798]) ...
    mean([4.1796    4.1794]) ...
    mean([4.1813    4.1813]) ...
    mean([4.1850    4.1850]) ...
    mean([4.1901    4.1901]) ...
    mean([4.1968    4.1968]) ...
    mean([4.2052    4.2052]) ...
    mean([4.2136    4.2136]) ...
    ]; % Specific heat at constant pressure table [kJ/(kg*K) -> J/(kg*K)]

k_TLU = 1E-3*[ ...
    mean([561.0400  561.0900]) ...
    mean([580.0200  580.0700]) ...
    mean([598.4400  598.4800]) ...
    mean([615.4800  615.5200]) ...
    mean([630.6000  630.6400]) ...
    mean([643.6100  643.6100]) ...
    mean([654.3900  654.3900]) ...
    mean([663.1300  663.1300]) ...
    mean([670.0200  670.0200]) ...
    mean([675.2700  675.2700]) ...
    mean([680.5200  680.5200]) ...
    ]; % Thermal conductivity table [mW/(m*K) -> W/(m*K)]

mu_TLU = 1E-3*[ ...
    mean([1.79134166000000	1.79104172000000]) ...
    mean([1.30570817000000	1.30550823000000]) ...
    mean([1.00139424000000	1.00139424000000]) ...
    mean([0.797077360000000	0.796977800000000]) ...
    mean([0.652569940000000	0.652569940000000]) ...
    mean([0.546364000000000	0.546364000000000]) ...
    mean([0.465938480000000	0.465938480000000]) ...
    mean([0.403538060000000	0.403538060000000]) ...
    mean([0.354026740000000	0.354026740000000]) ...
    mean([0.314108620000000	0.314108620000000]) ...
    mean([0.274792080000000	0.274792080000000]) ...
    ]; % Dynamic viscosity table [cP (Centipoise) -> Pa*S]

Pr_TLU = [
    mean([13.4736964762477	13.4686439490420]) ...
    mean([9.44467195482053	9.44151149977244]) ...
    mean([7.00176110960497	7.00079117123379]) ...
    mean([5.41343840991746	5.41202204386535]) ...
    mean([4.32521617701237	4.32473488398453]) ...
    mean([3.54952811982412	3.54952811982412]) ...
    mean([2.97980185944162	2.97980185944162]) ...
    mean([2.54982405441769	2.54982405441769]) ...
    mean([2.21751503303185	2.21751503303185]) ...
    mean([1.95609099889526	1.95609099889526]) ...
    mean([1.70143994046905	1.70143994046905]) ...
    ]; % Prandtl number table []


% outputs
rho = interp1(T_TLU, rho_TLU, T, 'makima', 'extrap');
cp  = interp1(T_TLU, cp_TLU, T, 'makima', 'extrap');
k   = interp1(T_TLU, k_TLU, T, 'makima' ,'extrap');
mu  = interp1(T_TLU, mu_TLU, T, 'makima' ,'extrap');
Pr  = interp1(T_TLU, Pr_TLU, T, 'makima' ,'extrap');
end