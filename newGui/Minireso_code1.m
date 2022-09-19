
close all;
clear all;

% ========================================
% Intitulé : Minireso_code1
% Auteur : Arnaud Fernandez
% Version du : 11/02/2020
%
% Fibre dopée : IXblue Rad [1 cm]
% Laser signal : [1530 - 1565 nm]
% Laser pompe : 976 nm 
% 
% Affiche le gain petit signal  
%
% ========================================

% Load parameters  iXBlueRAD ----------------------------------------------

root='./RadIX/iXBlueRADpart1';
[lmbd1,Y1]=datacq1(root);  % [nm,m2] lambda,Sap

figure(1)
subplot(111)
fs=20; % fontsize
plot(lmbd1,Y1*10^25,'b-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
axis([900 1050 0 2.5]);
legend('\sigma_{ap}');
grid on;
xlabel ('\lambda [nm]');
ylabel ('[10^{-25}.m^2]');
%title('IX fiber : Rad hard');

root='./RadIX/iXBlueRADpart2';
[lambda,S_as,S_es]=datacq2(root); % [nm,m2,m2] lambda,Sas,Ses
%lambda=lambda*1e-9; % (m)

figure(2)
subplot(111)
fs=20; % fontsize
plot(lambda,S_as*10^25,'b-',lambda,S_es*10^25,'r-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
axis([1450 1650 0  inf]);
legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('\lambda [nm]');
ylabel ('[10^{-25}.m^2]');
title('IX fiber : Rad hard');

N=1.4e25;	 % (ions.m^(-3))	Densité d’ions Er3+
% N2	% (ions.m(-3))	Densité d'ions sur le niveau excité
% Sap=2.3e-25; % (m2)	Section efficace d’absorption de la pompe   @ 980 nm
Sep=0;       % (m2)	Section efficace d’émission de la pompe @ 980 nm 
Gp=0.88;     % (sd) Facteur de recouvrement entre le mode de la pompe et l’aire dopée
Gs=0.52;	 % (sd) Facteur de recouvrement entre le mode du signal et l’aire dopée
Aeff=8.043e-12; % (m2)	Aire de la zone dopée
ap=1.3816e-3;% (m-1)	Pertes linéiques de la pompe
as=1.3816e-3;% (m-1)	Pertes linéiques de la sonde
Beta_sp=10^-6; % (sd) mode coupling factor

Tau=10e-3;   % (s) Metastable level lifetime
c=3e8;       % (m/s) light speed in vacuum 
ng=1.5;      % (sd) fiber group index
vg=c/ng;     % (m/s) group light speed 



% End load parameters  iXBlueRAD ----------------------------------------------

lbd_p=976*1e-9;                      % (m) lambda pump
lbd_s=1550*1e-9;                     % (m) lambda signal
w_s=2*pi*vg/lbd_s;                   % (rad.s-1) signal pulsation   
w_p=2*pi*vg/lbd_p;                   % (rad.s-1) pump pulsation
dw=(abs(lambda(1)-lambda(2))*(w_s)^2)/(c*2*pi);  % (rad.s-1) step

[Sap] = findSap(lbd_p,lmbd1*1e-9,Y1);% (m2)
[Sas,Ses]=findSasSes(lbd_s,lambda*1e-9,S_as,S_es);     % (m2)
Sap=2.254e-25;
Sas=2.6588e-25;
Ses=3.9142e-25;

L=0.01;                               % (m) FP resonator length 

Ra=99.996 ;                            % (%) mirrors power reflectivity 
Rb=99.996 ;                            % (%) mirrors power reflectivity 

h=6.62e-34;                           % (J.s) h Planck
hbar=h/(2*pi);
%Beta_sp=Gs;                          % (sd) fraction of spontaneous emission coupled into the mode
%r_sp=N2*??/Tau;                      % (s-1.m-3)/(rad.s-1) spectral density of spontaneous emission 
%Dsp=Beta_sp*hbar*w_s*r_sp;           % (W.m-2).m-1/(rad.s-1) spatial-spectral density of spontaneous emission intensity 

%=================================
% Gain petit signal de la fibre dopée
%=================================

Pp=1e-3; % (W) pump power
GdB= smallsig_gain(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N);

%GdB=0.0.691 dB @ Pp=1 mW;
%GdB=0.1908 dB @ Pp=5 mW;
%GdB=0.2212 dB @ Pp=15 mW;
%GdB=0.2362 dB @ Pp=150 mW;
%GdB=0.2375 dB @ Pp=500 mW;

Pp_min=0.000001e-3; % (W) min pump power
Pp_max=10e-3; % (W) max pump power
nbpt=1000;
Ppvec=linspace(Pp_min,Pp_max,nbpt);
Yvec=linspace(Pp_min,Pp_max,nbpt);
XdBvec=linspace(Pp_min,Pp_max,nbpt);
GdBvec=linspace(Pp_min,Pp_max,nbpt);
gthdBvec=linspace(Pp_min,Pp_max,nbpt);
Dnuvec=linspace(Pp_min,Pp_max,nbpt);
Qvec=linspace(Pp_min,Pp_max,nbpt);
Fvec=linspace(Pp_min,Pp_max,nbpt);
%Y=YvsPp(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N,Ra,Rb,hbar,w_s,Beta_sp,S_es,dw)

for ii=1:nbpt
    Pp=Ppvec(1,ii);
    Yvec(1,ii)=YvsPp(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N,Ra,Rb,hbar,w_s,Beta_sp,S_es,dw);
    [XdBvec(1,ii),gthdBvec(1,ii),Dnuvec(1,ii),Qvec(1,ii),Fvec(1,ii)] = XdBvsPp(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N,Ra,Rb,Yvec(1,ii),c,ng,w_s);
    GdBvec(1,ii)= smallsig_gain(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N);
end    

figure(3)
subplot(111)
fs=20; % fontsize
plot(Ppvec,Yvec,'b-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
%axis([1450 1650 0  inf]);
%legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('P pump [W]');
ylabel ('Iint');
title('IX fiber : Rad hard');

figure(4)
subplot(111)
fs=20; % fontsize
plot(Ppvec,XdBvec,'b-',Ppvec,gthdBvec,'r-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
%axis([1450 1650 0  inf]);
%legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('P pump [W]');
ylabel ('\Delta gain');
title('IX fiber : Rad hard');

figure(5)
subplot(111)
fs=20; % fontsize
plot(Ppvec,GdBvec,'b-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
%axis([1450 1650 0  inf]);
%legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('P pump [W]');
ylabel ('Small sig gain');
title('IX fiber : Rad hard');

figure(6)
subplot(111)
fs=20; % fontsize
plot(Ppvec,Qvec,'b-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
%axis([1450 1650 0  inf]);
%legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('P pump [W]');
ylabel ('Q factor');
title('IX fiber : Rad hard');

function [GdB] = smallsig_gain(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N)

nu_p=w_p/(2*pi);                    % 
Phi_Pth=1/(Tau*Sap);                % (s-1.m-2) pump photons flux at threshold
Phi_P=Pp*Gp/(h*nu_p*Aeff);          % (s-1.m-2) pump photons flux
Yp=Phi_P/Phi_Pth;                   % (sd) normalized pump photons flux
g0=Ses*N;                           % (m-1) maximum gain
g=g0*(Yp-1)/(Yp+1);                 % (m-1) small signal gain
GdB=10*log10(exp(g*L));             % (dB) small signal gain (single -pass)

end

function [XdB,gthdB,Dnu,Q,F] = XdBvsPp(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N,Ra,Rb,Y,c,ng,w_s)

nu_p=w_p/(2*pi);                    % 
Phi_Pth=1/(Tau*Sap);                % (s-1.m-2) pump photons flux at threshold
Phi_P=Pp*Gp/(h*nu_p*Aeff);          % (s-1.m-2) pump photons flux
Yp=Phi_P/Phi_Pth;                   % (sd) normalized pump photons flux
g0=Ses*N;                           % (m-1) maximum gain
g=g0*(Yp-1)/(Yp+1);                 % (m-1) small signal gain
%GdB=10*log10(exp(g*L));             % (dB) small signal gain (single -pass)
ra=sqrt(Ra/100);
rb=sqrt(Rb/100);                    % (sd) amplitude reflectivity
gth=(1/L)*log(1/(ra*rb));           % (m-1) threshold gain
Kth=gth/g0;
X=g0*L*(Kth-(Yp-1)/(Yp+1+2*Y));     % (sd) delta(g)L 
gthdB=10*log10(exp(gth*L));                  % (dB) threshold gain
XdB=10*log10(exp(X));               % (dB) delta gain gain seuil laser - gain intra-cavité 
Dnu=(c/(2*pi*ng*L))*X/sqrt(1-X);              % (Hz) FWHM of cavity resonance
FSR=c/(2*pi*ng*L);
Q=w_s/(2*pi*Dnu);                          % (sd) Q factor
F=FSR/Dnu;% (sd) Finesse
end

function [Y] = YvsPp(w_p,Tau,Sap,h,Gp,Aeff,Pp,L,Ses,N,Ra,Rb,hbar,w_s,Beta_sp,S_es,dw)

nu_p=w_p/(2*pi);                    % 
Phi_Pth=1/(Tau*Sap);                % (s-1.m-2) pump photons flux at threshold
Phi_P=Pp*Gp/(h*nu_p*Aeff);          % (s-1.m-2) pump photons flux
Yp=Phi_P/Phi_Pth;                   % (sd) normalized pump photons flux
g0=Ses*N;                           % (m-1) maximum gain
g=g0*(Yp-1)/(Yp+1);                 % (m-1) small signal gain
N2=0.5*((g/Ses)+N);                 % (m-3)
ra=sqrt(Ra/100);
rb=sqrt(Rb/100);                    % (sd) amplitude reflectivity
gth=(1/L)*log(1/(ra*rb));           % (m-1) threshold gain
Kth=gth/g0;
gw=Ses/(trapz(S_es)*dw);
r_sp=N2*gw/Tau;
%r_sp=1;
K=Beta_sp*hbar*w_s*r_sp*4*L*pi/Tau;
A=g0*L*Kth;
B=g0*L*(Kth*Yp+Kth-Yp-1-K);
C=-K*(Yp+1);
pol=[A B C];
sol=max(roots(pol));
Y=sol;

end



%--------------------------------------------------------------------------
% RK4
            % Runge-Kutta
            k1=dttc*(gL-g(:,2))-(dtEsat*g(:,2)).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
            g1=g(:,2)+0.5*k1;
            k2=dttc*(gL-g1)-(dtEsat*g1).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
            g2=g(:,2)+0.5*k2;
            k3=dttc*(gL-g2)-(dtEsat*g2).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
            g3=g(:,2)+k3;
            k4=dttc*(gL-g3)-(dtEsat*g3).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
                  
            g(:,1)=g(:,2)+(k1+2*k2+2*k3+k4)/6;  
 %--------------------------------------------------------------------------




% Ps_in_dBm=-40;                          % (dBm) input probe power
% Ps_in_W=(10.^(Ps_in_dBm/10))*(1e-3);    % (W) input probe power
% Pp_in_W=100e-3;                         % (W) input pump power 
% Pp_in=Pp_in_W;
% Lbd_sig=linspace(1525e-9,1565e-9,50);   % (m) input signal wavelength 
% 
% %=================================
% % Propagation
% %=================================
% 
% for LL=1:size(Lbd_sig,2)
%     
% lbd_s=Lbd_sig(1,LL);              % (m) input signal wavelength
% [Sas,Ses]=findSasSes(lbd_s,lambda*1e-9,S_as,S_es);     % (m2)
%     
% for kk=1:size(Lfib,2)
% L=Lfib(1,kk);        % (m) longueur fibre dopée 
% %=================================
% % Domaines temporel et fréquentiel
% %=================================
% 
% % SPATIAL (XX) AND SPATIAL FREQUENCY (FF) ARRAY PARAMETERS
% 
% dz=0.2;                 % (m) spatial step index 
% ndz=ceil(L/dz);         % (sd)
% Lz=linspace(1,L,ndz);
% 
% dt=dz/vg;               % (s) temporal step
% % tmax =(ntot-1)*dt/2;            % (ns)
% % dv = 1/(2*tmax);                % Step fréquentiel (GHz)
% % TT = [-ntot/2:(ntot/2)-1]*dt;   % Vecteur temps (s)
% % FF = [-ntot/2:(ntot/2)-1]*dv;   % Vecteur freq  (Hz)
% % tps=linspace(0,(size(TT,2)-1)*dt,size(TT,2));
% 
% tps_simul=0.25*Tau;  % (s)
% tps_fib=L/vg;      % (s)
% ntot=ceil(tps_simul/dt); % (s)c
% %pause
% 
% % Conditions initiales
% N2=zeros(2,ndz);
% N1=zeros(2,ndz);
% N1(1,:)=N;
% Pp=zeros(1,ndz);Ps=zeros(1,ndz);
% Pp_out=zeros(1,ntot);Ps_out=zeros(1,ntot);   
%     
% Ps_in=Ps_in_W;
% 
% 
% fprintf(1, '\nSimulation running...      ');
% 
% for jj=1:ntot
% 
% Pp(1,1)=Pp_in;
% Ps(1,1)=Ps_in;
% 
%     for ii=2:ndz
% 
%     % Equation de propagation de la pompe----------------------------------
%     B=Gp*(Sep*N2(1,ii-1)-Sap*N1(1,ii-1))-ap;
%     Pp(1,ii)=Pp(1,ii-1)*exp(B*dz); 
% 
%     % Equation de propagation de la sonde----------------------------------
%     C=Gs*(Ses*N2(1,ii-1)-Sas*N1(1,ii-1))-as;    %+2*Ses*N2(1,:)*(h*c^2)/(lbd_s^3 )*dlbd_s;
%     Ps(1,ii)=Ps(1,ii-1)*exp(C*dz); 
% 
%     end
% 
%     Ps_out(1,jj)=Ps(1,ndz); % output sonde
%     Pp_out(1,jj)=Pp(1,ndz); % output pump
% 
% 
%     % Equation d'évolution des porteurs------------------------------------
%     N2(2,:)=(-N2(1,:)/Tau+((Gp*lbd_p)/(h*c*A))*(Sap*N1(1,:)-Sep*N2(1,:)).*Pp(1,:)+((Gs*lbd_s)/(h*c*A))*(Sas*N1(1,:)-Ses*N2(1,:)).*Ps(1,:))*dt+N2(1,:);
%     
%     N2(1,:)=N2(2,:);
%     N1(1,:)=N-N2(1,:);
% 
%     fprintf(1, '\b\b\b\b\b\b%5.2f%%', jj * 100.0 /ntot );
%     
% end
% 
% %N2_final(kk,:)=N2(2,:);
% % Gain_s(LL,kk)=10*log10(Ps_out(1,size(Ps_out,2))/Ps_in);
% % Gain_p(LL,kk)=10*log10(Pp_out(1,size(Pp_out,2))/Pp_in);
% 
% Gain_s(LL,kk)=10*log10(Ps_out(1,size(Ps_out,2))/Ps_in);
% Gain_p(LL,kk)=10*log10(Pp_out(1,size(Pp_out,2))/Pp_in);
% 
% end
% 
% end
% 
% save('test15.mat');
% 
% 
% % % %=================================
% % % % Affichage
% % % %=================================
% % % 
% % % fprintf('  ++++++  SIGNAL ENTRANT: ++++++   \n');
% % % fprintf('Puissance moyenne sonde (dBm): %e\n',10*log10(Ps_in*1e3));
% % % fprintf('Puissance moyenne pompe (dBm): %e\n',10*log10(Pp_in*1e3));
% % % 
% % % Lz=linspace(1,L,ndz);
% % % figure(1)
% % % plot(Lz,N2(1,:),Lz,N1(1,:))
% % % legend('N2','N1')
% % % grid on
% % % 
% % % tps=linspace(0,ntot*dt,ntot);
% % % figure(2)
% % % plot(tps,Pp_out,tps,Ps_out,'LineWidth',1.5);
% % % set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% % % legend('Pump','Sonde')
% % % ylabel ('Power [W]');
% % % xlabel ('Time [s]');  
% % % grid on
% % % 
% % % figure(3)
% % % plot(tps,(Ps_out/Ps_in),'LineWidth',1.5);
% % % set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% % % legend('Sonde')
% % % ylabel ('Gain [dB]');
% % % xlabel ('Time [s]');  
% % % grid on