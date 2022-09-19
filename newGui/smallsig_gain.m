function [outputArg1,outputArg2] = smallsig_gain(inputArg1,inputArg2)


Phi_Pth=1/(Tau*Sig_P);              % (s-1.m-2) pump photons flux at threshold
Phi_P=Pp*Gam_P/(h*nu_p*Aeff);       % (s-1.m-2) pump photons flux
Yp=Phi_P/Phi_Pth;                   % (sd) normalized pump photons flux
g0=Sig_s*N;                         % (m-1) maximum gain
g=g0*(Yp-1)/(Yp+1);                 % (m-1) small signal gain
GdB=10*log10(exp(g*L));             % (dB) small signal gain (single -pass)


end

