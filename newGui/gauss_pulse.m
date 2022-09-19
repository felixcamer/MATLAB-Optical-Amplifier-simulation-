% *****************************************************************
% Generation of a gaussian pulse 
% Author : Arnaud Fernandez
% Last modification : 24/10/2014
% *****************************************************************

function [u0] = gauss_pulse(m_Gauss,C,PeakPower,Tp,t)


T0=Tp/(2*sqrt(log(2)));         % (ps) 1/e=0.36 halfwidth
ni0=1/T0;
u0=sqrt(PeakPower)*exp(-0.5*(1+1i*C)*((t*ni0).^(2.0*m_Gauss)));  % (W^0.5) Envelop

end

