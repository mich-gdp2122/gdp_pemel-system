function y = blend(y1,y2,x1,x2,x)
% Blend between two values y1 and y2 as x varies between x1 and x2. The
% blending function ensures y is continuous and differentiable in the
% range x1 <= x <= x2.
% 
% Copyright 2017-2018 The MathWorks, Inc.
% 
% [Source: simscape.function.blend.scc]
%
%
%%
u = (x-x1)/(x2-x1);
t = 3*u^2 - 2*u^3;
y = (1-t)*y1 + t*y2;
end