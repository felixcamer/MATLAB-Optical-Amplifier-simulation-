
close all;
clear all;

% =========================================================================
% Name : dopedfiber_TelecomIX_RadIX
% Author : Arnaud Fernandez
% Version du : 20/05/2020
%
% Doped fiber : IXblue Rad ou Telecom[qqs m]
% Signal laser : [1530 - 1565 nm]
% Pump laser : around 976 nm 
% 
% Calculates : adapted to multiple Pump power and fixed fiber length and
% Signal laser power
% forward and backward ASE ; cw signal laser amplification ; NF 
% 
% Physical properties integrated by the model :
% cw pump depletion along the doped fiber ; Er3+ ions density dependance over length ;
% Er3+ transition in glass is homogeneously broadened.
%
% =========================================================================


flag_Sae=1;         % 1 : display of fiber Sae parameters ; 0: don't display
flag_RTevo=1;       % 1 : display power field evolution in the fiber-length ; 0: don't display
ASEcalc=1;          % 1 : ASE included in the model ; 0 : not included
choix_fibre=0;      % 0 : Telecom fiber ; 1 : Rad fiber  

switch choix_fibre
    
    case 0
        
% Load parameters  iXBlueTelecom ----------------------------------------------

root='./TelecomIX/';

[num,txt,raw] = xlsread( fullfile(root, 'iXBlueTelpart2.xlsx') );

[lmbd1,Y1]=datacq3(root);  % [nm,m2] lambda,Sap

if flag_Sae==1
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
end

[lambda,S_as,S_es]=datacq4(root); % [nm,m2,m2] lambda,Sas,Ses
lambda=lambda*1e-9; % (m)

if flag_Sae==1
figure(2)
subplot(111)
fs=20; % fontsize
plot(lambda*1e9,S_as*10^25,'b-',lambda*1e9,S_es*10^25,'r-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
axis([1450 1650 0  inf]);
legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('\lambda [nm]');
ylabel ('[10^{-25}.m^2]');
title('IX fiber : Telecom');
end

N=5e24;	     % (ions.m^(-3))	Densité d’ions Er3+
Sep=0;       % (m2)	Section efficace d’émission de la pompe @ 980 nm 
Gp=0.87;     % (sd) Facteur de recouvrement entre le mode de la pompe et l’aire dopée
Gs=0.52;	 % (sd) Facteur de recouvrement entre le mode du signal et l’aire dopée
Aeff=5.73e-12;  % (m2)	Aire de la zone dopée
ap=1.612e-3; % (m-1)	Pertes linéiques de la pompe
as=1.612e-3; % (m-1)	Pertes linéiques de la sonde
Tau=10e-3;   % (s) Durée de vie du niveau métastable
ng=1.5;      % (sd) fiber group index

% End load parameters  iXTelecom ------------------------------------------
 
    case 1

% Load parameters  iXBlueRAD ----------------------------------------------

root='./RadIX/iXBlueRADpart1';
[lmbd1,Y1]=datacq1(root);  % [nm,m2] lambda,Sap

if flag_Sae==1
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
end

root='./RadIX/iXBlueRADpart2';
[lambda,S_as,S_es]=datacq2(root); % [nm,m2,m2] lambda,Sas,Ses
lambda=lambda*1e-9; % (m)

if flag_Sae==1
figure(2)
subplot(111)
fs=20; % fontsize
plot(lambda*1e9,S_as*10^25,'b-',lambda*1e9,S_es*10^25,'r-','LineWidth',1.5); 
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
axis([1450 1650 0  inf]);
legend('\sigma_{a}','\sigma_{e}');
grid on;
xlabel ('\lambda [nm]');
ylabel ('[10^{-25}.m^2]');
title('IX fiber : Rad hard');
end

N=1.4e25;       % (ions.m^(-3))	Densité d’ions Er3+
Sep=0;          % (m2)	Section efficace d’émission de la pompe @ 980 nm 
Gp=0.88;        % (sd) Facteur de recouvrement entre le mode de la pompe et l’aire dopée
Gs=0.52;        % (sd) Facteur de recouvrement entre le mode du signal et l’aire dopée
Aeff=8.043e-12; % (m2)	Aire de la zone dopée
ap=1.3816e-3;   % (m-1)	Pertes linéiques de la pompe
as=1.3816e-3;   % (m-1)	Pertes linéiques de la sonde
Tau=10e-3;      % (s) Metastable level lifetime
ng=1.5;         % (sd) fiber group index

end

c=3e8;          % (m/s) light speed in vacuum 
vg=c/ng;        % (m/s) group light speed
h=6.62e-34;     % (J.s) h Planck
hbar=h/(2*pi);

pause

%==========================================================================
% Initial conditions

Lf=20;                                 % (m) doped fiber length
lbd_p=976*1e-9;                        % (m) lambda pump
nu_p=c/lbd_p;                          % (Hz) pump carrier frequency 
[Sap] = findSap(lbd_p,lmbd1*1e-9,Y1);  % (m2) 
Ps_in_dBm=-40;                         % (dBm) input signal power
Ps_in_W=(10.^(Ps_in_dBm/10))*(1e-3);   % (W) input signal power
% Pp_in_W=80e-3;                        % (W) input pump power 
lbd_s=1550e-9;                           % (m) input signal wavelength 
nu_s=c/lbd_s;                            % (Hz) signal carrier frequency 
dlbd=abs(lambda(2)-lambda(1));           % delta lambda (m)
dnu=(c/(lbd_s^2))*dlbd;                  % (Hz) delta nu
nu=c./lambda;                            % (Hz) freq. vector
[Sas,Ses]=findSasSes(lbd_s,lambda,S_as,S_es);     % (m2)
Cas=(Tau*Sas*Gs)/(Aeff*h*nu_s);
Cap=(Tau*Sap*Gp)/(Aeff*h*nu_p); 
C=(Tau*Gs)/(Aeff*h);
vCa=(lambda'/c).*(S_as)';
vCae=(lambda'/c).*(S_as+S_es)';
Ces=(Tau*Ses*Gs)/(Aeff*h*nu_s);
Cep=(Tau*Sep*Gp)/(Aeff*h*nu_p);         
dz=0.001;                % (m) spatial step index 
ndz=ceil(Lf/dz);         % (sd) nb of spatial steps
Lz=linspace(0,Lf,ndz);   % (m) spatial vector describing the fiber
dt=dz/vg;                % (s) temporal step
nlbd=size(lambda,2);     % (sd) nb of spectral slices   
tps_fib=Lf/vg;           % (s) fiber propagation time
tps_simul=200*tps_fib;     % (s) simulation time 
ntot=ceil(tps_simul/dt); % (sd) nb of temporal steps
time=linspace(0,tps_simul,ntot); % (s) time simulation vector


%PpL=zeros(1,ntot);PsL=zeros(1,ntot);Am0=zeros(1,ntot);ApL=zeros(1,ntot); 
nbsim=50;
tabPp_in=linspace(0.001,0.150,nbsim); % (W) pump power array
%sim_nb=linspace(8,27,20);          % sim number (small signal gain with ASE)
sim_nb=linspace(68,117,nbsim);         % sim number (small signal gain without ASE)

for kk=1:size(tabPp_in,2)          % (W) input pump power 
Pp_in_W=tabPp_in(1,kk);

% Vector declaration and initialization -----------------------------------

N2=zeros(1,ndz);            % (m-3) excitated state
N1=zeros(1,ndz);N1(1,:)=N;  % (m-3) fundamental state
Pp=zeros(1,ndz); % (W) pump power
Ps=zeros(1,ndz); % (W) signal power
Am=zeros(nlbd,ndz); % (W) backward ASE
Ap=zeros(nlbd,ndz); % (W) forward ASE
% Outgoing power : PpL : pump @ z=Lf  (W) ; PsL : signal @ z=Lf  (W); 
% Am0 : backward ASE @ z=0  (W); ApL : forward ASE @ z=Lf  (W)


%==========================================================================
% Propagation equations

fprintf(1, '\nSimulation running...      ');
tstart = tic; % start time watch
% Regroupment de constantes utiles
      SeasV=S_es'+S_as';
      NSasV=N*S_as';
      ASEct=(S_es')*Gs*h.*(nu')*dnu;

% Ap1=zeros(nlbd,ndz);Ap2=zeros(nlbd,ndz);Ap3=zeros(nlbd,ndz);
% Am1=zeros(nlbd,ndz);Am2=zeros(nlbd,ndz);Am3=zeros(nlbd,ndz);
% k1p=zeros(nlbd,ndz);k2p=zeros(nlbd,ndz);k3p=zeros(nlbd,ndz);
% k1m=zeros(nlbd,ndz);k2m=zeros(nlbd,ndz);k3m=zeros(nlbd,ndz);
k1=zeros(nlbd,ndz);k2=zeros(nlbd,ndz);k3=zeros(nlbd,ndz);

%Pp=(gauss_pulse(10,0,Pp_in_W,0.99*tps_simul,time)).^2; % (W^0.5)
Trise=tps_fib;      % (s) rise time for a slow turn-on of input signal and pump laser 
PpIN=Pp_in_W*(1-exp(-time/Trise));
PsIN=Ps_in_W*(1-exp(-time/Trise));

count=0;jj=0;
%for jj=1:ntot
while count < 100
    
    jj=jj+1;
    % Boundary conditions
      Ap(:,1)=0; Am(:,ndz)=0;
      Pp(1,1)=PpIN(1,jj); Ps(1,1)=PsIN(1,jj);
    
    % Population equation ------------------------------------------------- 
      N2(1,:)=N*(Cas*Ps(1,:)+Cap*Pp(1,:)+C*sum(vCa.*(Ap+Am),1))./(1+(Cas+Ces)*Ps(1,:)+(Cap+Cep)*Pp(1,:)+C*sum(vCae.*(Ap+Am),1));  
    
    % Waves propagation ---------------------------------------------------
      
      Pp(1,2:ndz)=Pp(1,1:ndz-1).*exp(((N2(1,1:ndz-1)*(Sep+Sap)-N*Sap)*Gp-ap)*dz); 
      Ps(1,2:ndz)=Ps(1,1:ndz-1).*exp(((N2(1,1:ndz-1)*(Ses+Sas)-N*Sas)*Gs-as)*dz);
      
% %      k1=dz*(N2(1,1:ndz-1).*(S_es'+S_as')-N*S_as')*Gs.*Ap(:,1:ndz-1)+N2(1,1:ndz-1).*(S_es')*Gs*h.*(nu')*dnu-as*Ap(:,1:ndz-1);
%        k1p=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap(:,1:ndz-1);
%        Ap1=Ap(:,1:ndz-1)+0.5*k1p;
%        k2p=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap1(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap1(:,1:ndz-1);
%        Ap2=Ap(:,1:ndz-1)+0.5*k2p;
%        k3p=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap2(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap2(:,1:ndz-1);
%        Ap3=Ap(:,1:ndz-1)+k3p;
%        k4p=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap3(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap3(:,1:ndz-1);
%        Ap(:,2:ndz)=Ap(:,1:ndz-1)+(1/6)*(k1p+2*k2p+2*k3p+k4p);
%               
%        k1m(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am(:,2:ndz);
%        Am1(:,2:ndz)=Am(:,2:ndz)+0.5*k1m(:,2:ndz);
%        k2m(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am1(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am1(:,2:ndz);
%        Am2(:,2:ndz)=Am(:,2:ndz)+0.5*k2m(:,2:ndz);
%        k3m(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am2(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am2(:,2:ndz);
%        Am3(:,2:ndz)=Am(:,2:ndz)+k3m(1,2:ndz);
%        k4m(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am3(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am3(:,2:ndz);
%        Am(:,1:ndz-1)=Am(:,2:ndz)+(1/6)*(k1m(:,2:ndz)+2*k2m(:,2:ndz)+2*k3m(:,2:ndz)+k4m(:,2:ndz));     

if ASEcalc==1
              
       k1=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap(:,1:ndz-1);
       %Ap1=Ap(:,1:ndz-1)+0.5*k1p;
       k2=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+0.5*k1)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+0.5*k1);
       %Ap2=Ap(:,1:ndz-1)+0.5*k2p;
       k3=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+0.5*k2)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+0.5*k2);
       %Ap3=Ap(:,1:ndz-1)+k3p;
       k4=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+k3)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+k3);
       Ap(:,2:ndz)=Ap(:,1:ndz-1)+(1/6)*(k1+2*k2+2*k3+k4);
              
       k1(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am(:,2:ndz);
       %Am1(:,2:ndz)=Am(:,2:ndz)+0.5*k1m(:,2:ndz);
       k2(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+0.5*k1(:,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+0.5*k1(:,2:ndz));
       %Am2(:,2:ndz)=Am(:,2:ndz)+0.5*k2m(:,2:ndz);
       k3(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+0.5*k2(:,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+0.5*k2(:,2:ndz));
       %Am3(:,2:ndz)=Am(:,2:ndz)+k3m(1,2:ndz);
       k4(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+k3(1,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+k3(1,2:ndz));
       Am(:,1:ndz-1)=Am(:,2:ndz)+(1/6)*(k1(:,2:ndz)+2*k2(:,2:ndz)+2*k3(:,2:ndz)+k4(:,2:ndz)); 
else
       Ap(:,2:ndz)=0;
       Am(:,1:ndz-1)=0;
end             
%             % Runge-Kutta
%             k1=dttc*(gL-g(:,2))-(dtEsat*g(:,2)).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
%             g1=g(:,2)+0.5*k1;
%             k2=dttc*(gL-g1)-(dtEsat*g1).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
%             g2=g(:,2)+0.5*k2;
%             k3=dttc*(gL-g2)-(dtEsat*g2).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
%             g3=g(:,2)+k3;
%             k4=dttc*(gL-g3)-(dtEsat*g3).*((abs(Ep(:,2))).^2+(abs(Ec(:,2))).^2);
%                   
%             g(:,1)=g(:,2)+(k1+2*k2+2*k3+k4)/6;
            
       
    % Output waves -------------------------------------------------------- 
      PpL(1,jj)=Pp(1,ndz);  % (W) outgoing pump 
      PsL(1,jj)=Ps(1,ndz);  % (W) outgoing signal
      if jj>ceil(Trise/dt)
         if PsL(1,jj)==PsL(1,jj-1)
             count=count+1;
         end    
      end    
      
      ApL(1,jj)=sum(Ap(:,ndz),1);  % (W) outgoing forward ASE @ z=L
      Am0(1,jj)=sum(Am(:,1),1);    % (W) outgoing backward ASE @ z=0
      
    % Realtime display ----------------------------------------------------
    if flag_RTevo==1
        
        fs=12;
%         figure(2)
%         subplot(211)
%         Color='b';
%         fs=12;
%         plot(Lz,10*log10(sum(Ap,1)*1e3),Lz,10*log10(sum(Am,1)*1e3),'LineWidth',1.5);
%         set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
%         grid on;
%         xlabel ('length (m)');
%         ylabel ('ASE power (dBm)');
%         legend('forward','backward')
%         
%         subplot(212)
%         Color='b';
%         Color2='r';
%         [AZ,H1,H2]=plotyy(Lz,10*log10(Ps*1e3),Lz,10*log10(Pp*1e3));
%         set(H1,'Color',Color,'LineWidth',1.5);
%         set(H2,'Color',Color2,'LineWidth',1.5);
%         set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
%         set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
%         % xlim(AZ(1),[-51,51])
%         % xlim(AZ(2),[-51,51])
% 
%         % ylim([0 1.01])
%         % set(AX,'Fontsize',fs,'FontName','times', ...
%         % 'xlim',[-limX limX],'xtick',(-20*stepX:stepX:20*stepX)); % Axe absisses
% 
%         set(AZ(1),'box','off')
%         xlabel(AZ(1),'length (m)','Fontsize',fs);
%         ylabel(AZ(1),'Signal (dBm)','Fontsize',fs);
%         set(AZ(2),'box','off')
%         ylabel(AZ(2),'Pump (dBm)','Fontsize',fs);
%         grid on;
        
        figure(3)
        subplot(111)
        Color='b';
        Color2='r';
        fs=12;
        plot(PsL*1e3,Color,'LineWidth',1.5);
        %plot(lambda*1e9,10*log10(Am(:,1)*1e3),Color2,'LineWidth',1.5);
        set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
        grid on;
        xlabel ('time (ns)');
        ylabel ('power (mW)');
        %legend('forward','backward')


    end    
    %----------------------------
         
      
      
      fprintf(1, '\b\b\b\b\b\b%5.2f%%', jj * 100.0 /ntot );
    
      end

fprintf('\n');
telapsed = toc(tstart);
fprintf('Elapsed rime (s) : %f\n',telapsed);       % (s) elapsed time

% Gain_s(1,:)=10*log10((PsL*1e3)/(Ps_in_W*1e3));
% Gain_p(1,:)=10*log10((PpL*1e3)/(Pp_in_W*1e3));

tmp='';
tmp=strcat('test',int2str(sim_nb(1,kk)));
save(tmp);
clear ApL Am0 PpL PsL

end


%time=linspace(1,tps_simul,ntot); % (s) time simulation vector

% figure(2)
% fs=12;
% subplot(111)
% Color='b';
% Color2='r';
% [AZ,H1,H2]=plotyy(time*1e9,Gain_s,time*1e9,Gain_p);
% set(H1,'Color',Color,'LineWidth',1.5);
% set(H2,'Color',Color2,'LineWidth',1.5);
% set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
% set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
% % xlim(AZ(1),[-51,51])
% % xlim(AZ(2),[-51,51])
% 
% % ylim([0 1.01])
% % set(AX,'Fontsize',fs,'FontName','times', ...
% % 'xlim',[-limX limX],'xtick',(-20*stepX:stepX:20*stepX)); % Axe absisses
% 
% set(AZ(1),'box','off')
% xlabel(AZ(1),'time (ns)','Fontsize',fs);
% ylabel(AZ(1),'Signal gain (dB)','Fontsize',fs);
% set(AZ(2),'box','off')
% ylabel(AZ(2),'Pump gain (dB)','Fontsize',fs);
% grid on;
% 
%         
% fs=12;
% figure(3)
% subplot(211)
% Color='b';
% Color2='r';
% fs=12;
% plot(Lz,10*log10(sum(Ap,1)*1e3),Color,Lz,10*log10(sum(Am,1)*1e3),Color2,'LineWidth',1.5);
% set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% grid on;
% xlabel ('length (m)');
% ylabel ('ASE power (dBm)');
% legend('forward','backward')
% 
% subplot(212)
% Color='b';
% Color2='r';
% [AZ,H1,H2]=plotyy(Lz,10*log10(Ps*1e3),Lz,10*log10(Pp*1e3));
% set(H1,'Color',Color,'LineWidth',1.5);
% set(H2,'Color',Color2,'LineWidth',1.5);
% set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
% set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
% % xlim(AZ(1),[-51,51])
% % xlim(AZ(2),[-51,51])
% 
% % ylim([0 1.01])
% % set(AX,'Fontsize',fs,'FontName','times', ...
% % 'xlim',[-limX limX],'xtick',(-20*stepX:stepX:20*stepX)); % Axe absisses
% 
% set(AZ(1),'box','off')
% xlabel(AZ(1),'length (m)','Fontsize',fs);
% ylabel(AZ(1),'Signal (dBm)','Fontsize',fs);
% set(AZ(2),'box','off')
% ylabel(AZ(2),'Pump (dBm)','Fontsize',fs);
% grid on;


%save('test00.mat'); % pump 300 mW


