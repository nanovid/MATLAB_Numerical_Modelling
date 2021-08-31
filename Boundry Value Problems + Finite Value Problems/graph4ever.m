function [x_out, u_out] = graph4ever(n)

x_out = linspace(0, n)

u_out = ((2-exp(1))./((exp(1) - 1).*exp(1))).*exp(2.*x_out) + ((exp(2) - 2)./((exp(1) - 1).*exp(1))).*exp(x_out)
plot(x_out, u_out, '.r')
hold on
title("sshakim")
