%% Fourier Series
%
%   Written by Ted Kwan for Math 226B.
%
%   This function calculates the Fourier
%   series for -\Delta u=1 on the unit square
%   out to 110 coefficients
function [uj] = fourierPoisson(x,y)
    uj=0;
    for k=1:110
        %%% Calculate current summand
        %
        k21=(2*k-1);
        a=(8/((pi)^3));
        b=(sinh(k21*pi*(1-y))+sinh(k21*pi*y))./(sinh(k21*pi));
        c=sin((pi*(2*k-1))*x)./(k21^3);
        %%% Sum up coefficients
        %
        uj=uj+(a*b*c);
    end
    uj=((x*(1-x)-uj))/2;
end