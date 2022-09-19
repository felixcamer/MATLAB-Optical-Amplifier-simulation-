% =========================================================================
% Name : Propag_EDFA
% Author : Barkoum Betra Felix
% Version du : 03/06/2022
%
% Doped fiber : IXblue Rad ou Telecom [qqs m]
% Signal laser : [1530 - 1565 nm]
% Pump laser : around 976 nm 
% 
% Computes EDF propagation equation and rate equation
% eqs (6.3) to (6.7) p 157 ; EDFA [Becker]
%
% =========================================================================


% Vector declaration and initialization -----------------------------------

ndz=ceil(Lf/dz);         % (sd) nb of spatial steps
Lz=linspace(0,Lf,ndz);   % (m) spatial vector describing the fiber
dt=dz/vg;                % (s) temporal step
nlbd=size(lambda,2);     % (sd) nb of spectral slices    
tps_fib=Lf/vg;           % (s) fiber propagation time
tps_simul=200*tps_fib;     % (s) simulation time 
ntot=ceil(tps_simul/dt); % (sd) nb of temporal steps
time=linspace(0,tps_simul,ntot); % (s) time simulation vector

N2=zeros(ntot,ndz);            % (m-3) excitated state
N1=zeros(1,ndz);N1(1,:)=N;  % (m-3) fundamental state
Pp=zeros(1,ndz);    % (W) pump power
[ind_lbd]=findech(lambda,lbd_s);
Ps=zeros(1,ndz);    % (W) signal power
Am=zeros(nlbd,ndz); % (W) backward ASE
Ap=zeros(nlbd,ndz); % (W) forward ASE

counter=100;        % (nd) nb of required consecutive occurence of below tolerance matching
tth=40e-9;  % (s) inital threshold time for evaluation
tth_step = 40e-9; % (s) delta threshold time for evaluation
delT=10e-9; % (s) delta t
delinc=ceil(delT/dt); % (sd) delta samples
tol=5e-10;  % (W/ns) tolerance
% inc_prev=0;         % (nd) previous increment
% tol=0.0005;         % (%) relative error tolerance
% counter=200;        % (nd) nb of required consecutive occurence of below tolerance matching

PsL = zeros(1,ntot);    % (W) output signal power
ApL = zeros(1,ntot);
PpL = zeros(1,ntot);    % (W) output pump power
ApLs = zeros(1,ntot);
Am0 = zeros(1,ntot);

%==========================================================================
% Propagation equations

%fprintf(1, '\nSimulation running...      ');
tstart = tic; % start time watch
% Regroupment de constantes utiles
      SeasV=S_es'+S_as';
      NSasV=N*S_as';
      ASEct=(S_es')*Gs*h.*(nu')*dnu;

k1=zeros(nlbd,ndz);
k2=zeros(nlbd,ndz);
k3=zeros(nlbd,ndz);
k4=zeros(nlbd,ndz);

Trise=5e-9;      % (s) rise time for a slow turn-on of input signal and pump laser 
%Trise = 0;      % (s) rise time for a slow turn-on of input signal and pump laser 

i = 1;
Pp_in_W = PpINtab(i);
PpIN=Pp_in_W*(1-exp(-time/Trise));
PsIN=Ps_in_W*(1-exp(-time/Trise));
%PpIN=Pp_in_W*ones(1,ntot);
%PsIN=Ps_in_W*ones(1,ntot);

count=0; jj=0;
while count < counter
    
    if mod(jj,500) == 0
        disp(jj)
    end
    
    jj=jj+1;
    %tt(1,jj)=jj*dt;
    % Boundary conditions
      Ap(:,1)=0;
      Am(:,ndz)=0;
      Pp(1,1)=PpIN(1,jj);
      Ps(1,1)=PsIN(1,jj);
      N2(jj+1,:)= N2(1,:);
      
    % Population equation ------------------------------------------------- 
      N2(1,:)=N*(Cas*Ps(1,:)+Cap*Pp(1,:)+C*sum(vCa.*(Ap+Am),1))./(1+(Cas+Ces)*Ps(1,:)+(Cap+Cep)*Pp(1,:)+C*sum(vCae.*(Ap+Am),1));  
    
    % Waves propagation ---------------------------------------------------
      
      Pp(1,2:ndz)=Pp(1,1:ndz-1).*exp(((N2(1,1:ndz-1)*(Sep+Sap)-N*Sap)*Gp-ap)*dz); 
      Ps(1,2:ndz)=Ps(1,1:ndz-1).*exp(((N2(1,1:ndz-1)*(Ses+Sas)-N*Sas)*Gs-as)*dz);
      
    if ASEcalc==1
           % Runge-Kutta       
           k1=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*Ap(:,1:ndz-1)+N2(1,1:ndz-1).*ASEct-as*Ap(:,1:ndz-1);
           k2=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+0.5*k1)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+0.5*k1);
           k3=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+0.5*k2)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+0.5*k2);
           k4=dz*(N2(1,1:ndz-1).*SeasV-NSasV)*Gs.*(Ap(:,1:ndz-1)+k3)+N2(1,1:ndz-1).*ASEct-as*(Ap(:,1:ndz-1)+k3);
           Ap(:,2:ndz)=Ap(:,1:ndz-1)+(1/6)*(k1+2*k2+2*k3+k4);

           k1(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*Am(:,2:ndz)+N2(1,2:ndz).*ASEct-as*Am(:,2:ndz);
           k2(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+0.5*k1(:,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+0.5*k1(:,2:ndz));
           k3(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+0.5*k2(:,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+0.5*k2(:,2:ndz));
           k4(:,2:ndz)=dz*(N2(1,2:ndz).*SeasV-NSasV)*Gs.*(Am(:,2:ndz)+k3(1,2:ndz))+N2(1,2:ndz).*ASEct-as*(Am(:,2:ndz)+k3(1,2:ndz));
           Am(:,1:ndz-1)=Am(:,2:ndz)+(1/6)*(k1(:,2:ndz)+2*k2(:,2:ndz)+2*k3(:,2:ndz)+k4(:,2:ndz)); 
    else
           Ap(:,2:ndz)=0;
           Am(:,1:ndz-1)=0;
    end             
       
    % Output waves -------------------------------------------------------- 
      
      PpL(1,jj)=Pp(1,ndz);         % (W) outgoing pump 
      Ps0(1,jj)=Ps(1,1);           % (W) ingoing signal
      PsL(1,jj)=Ps(1,ndz);         % (W) outgoing signal
      ApL(1,jj)=sum(Ap(:,ndz),1);  % (W) outgoing forward ASE @ z=L
      ApLs(1,jj)=Ap(ind_lbd,ndz);  % (W) outgoing forward ASE @ z=L and lambda signal
      Am0(1,jj)=sum(Am(:,1),1);    % (W) outgoing backward ASE @ z=0
	  
	  
	  %-----------------------MINE-----------------------------------------
		Pp0(1,jj) = Pp(1,1);           % (W) ingoing signal
        Pp_evol(jj,:) = Pp(1,:);       % evolution of pump power
        Ps0(1,jj) = Ps(1,1);           % (W) ingoing signal
        Ps_evol(jj,:) = Ps(1,:);       % evolution of signal power 
        N2_evol(jj,:) = N2(1,:);       % evolution of N2 
		Ps_evol(jj,:) = Ps(1,:);       % evolution of signal power 
    
		N2_evol(jj,:) = N2(1,:);       % evolution of N2 
		
		Ap0(1,jj) = sum(Ap(:,1),1);    % (W) outgoing backward ASE @ z=0 (doit toujours etre = 0)
		ApLs(1,jj) = Ap(ind_lbd,ndz);% (W) outgoing forward ASE @ z=L and lambda signal
		ApL(1,jj) = sum(Ap(:,ndz),1);  % (W) outgoing forward ASE @ z=L
		ApL_lbd = Ap(:,ndz);           % (W) forward ASE (lbd) at z=L
		
		Ap_s_evol(jj,:) = Ap(ind_lbd,:); % evolution of ASE forward for signal wavelength
		
		Am0(1,jj) = sum(Am(:,1),1);      % (W) outgoing backward ASE @ z=0
		AmLs(1,jj) = Am(ind_lbd,ndz);  % (W) outgoing backward ASE @ z=0 and lambda signal (doit toujours etre = 0)
		AmL(1,jj) = sum(Am(:,ndz),1);    % (W) outgoing forward ASE @ z=L (doit toujours être = 0)
		Am0_lbd = Am(:,1);               % (W) backward ASE (lbd) at z=0
		
		Am_s_evol(jj,:) = Am(ind_lbd,:); % evolution of ASE backward for signal wavelength
		
		Ap_evol(jj,:) = sum(Ap,1);       % (W) Total Forward ASE (z)  
		Am_evol(jj,:) = sum(Am,1);       % (W) Total Backward ASE (z)  

		Gain_s = 10*log10((Ps(1,ndz)*1e3)/(Ps(1,1)*1e3));   % signal gain (dB)
		Gain_p = 10*log10((Pp(1,ndz)*1e3)/(Pp(1,1)*1e3));   % pump gain (dB)
		Gain_evol(jj,:) = 10*log10((Ps(1,:)*1e3)/(Ps(1,1)*1e3)); % signal gain evolution
		Gain_P_evol(jj,:) = 10*log10((Pp(1,:)*1e3)/(Pp(1,1)*1e3)); % signal gain evolution
		
		N2_evol_graph = N2_evol(jj,:)/N;    % N2/N normalized upper-state population
		  
	  %------------------------END MINE ------------------------------------------
      % Simulation ending evaluation -------------------------------------- 
      if time(1,jj)>tth % s   
        dPdt=abs(PsL(1,jj)-PsL(1,jj-delinc))/(abs(time(1,jj)-time(1,jj-delinc))*1e9);  % W/ns
          
            if  dPdt<tol
                count=count+1;
            end
          
            if  count == counter && (i + 1) <= size(PpINtab,2)
                
                Gain_s = 10*log10((Ps(1,ndz)*1e3)/(Ps(1,1)*1e3)); % signal gain (dB)
                Gain_p = 10*log10((Pp(1,ndz)*1e3)/(Pp(1,1)*1e3)); % pump gain (dB)
                telapsed = toc(tstart);
                Ps_in_W=(10.^(Ps_in_dBm/10))*(1e-3);    % (W) input signal power
                Pp_in_dBm=10*log10(Pp_in_W*1e3);  
                EDFA_RESTAB;
                tstart = tic;   %restart time
                
                i = i + 1;
                Pp_in_W = PpINtab(i);
                PpIN(1,jj+1:ntot) = PpINtab(i-1) + (PpINtab(i) - PpINtab(i-1))*(1-exp(-(time(1,jj+1:ntot)-time(1,jj+1))/Trise)); % increasing of previous value of PpIN 
                count = 0; %condition to stay inside loop while
                tth = tth_step + time(1,jj); % modification of tth for new simulation ending evaluation entry
                
                %verification of right simulation behaviour 
                % fs=12;
                % Color='b';
                % Color2='r';
                % figure(4)
                % yyaxis left
                % plot(time(1,1:jj)*1e9,PpIN(1,1:jj)*1e3,Color,'markerfacecolor',[0 0 1],'LineWidth',2.5,'LineStyle','-')
                % ylabel ('Input pump power (mW)');
                % yyaxis right
                % plot(time(1,1:jj)*1e9,PsL(1,1:jj),Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
                % set(gca,'Fontsize',fs,'FontName','times','XColor','k');
                % grid on;
                % xlabel ('time (ns)');
                % ylabel ('Output signal power (mW)');

            end
      end  
      % -------------------------------------------------------------------

    % Realtime display ----------------------------------------------------
    if flag_RTevo==1

        fs=12;
        Color='b';
        Color2='r';
        fs=12;
        figure(3)
        subplot(111)
        yyaxis left
        plot(tt*1e9,Ps0*1e3,Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
        ylabel ('input power (mW)');
        yyaxis right
        plot(tt*1e9,PsL*1e3,Color1,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
        set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
        grid on;
        xlabel ('time (ns)');
        ylabel ('output   power (mW)');

    end    
    
end

fprintf('\n');
telapsed = toc(tstart);
fprintf('Elapsed rime (s) : %f\n',telapsed);       % (s) elapsed time

Gain_s = 10*log10((Ps(1,ndz)*1e3)/(Ps(1,1)*1e3)); % signal gain (dB)
Gain_p = 10*log10((Pp(1,ndz)*1e3)/(Pp(1,1)*1e3)); % pump gain (dB)

disp('Gain sig = ');
disp(Gain_s);

disp('Gain pump = ');
disp(Gain_p);
%=====================================================================

% Gain_s_tt(1,:)=10*log10((PsL*1e3)/(Ps_in_W*1e3));
% Gain_p_tt(1,:)=10*log10((PpL*1e3)/(Pp_in_W*1e3));

% fs=12;
% Color='b';
% Color2='r';
% figure(3)
% subplot(211)
% plot(PpL,Gain_p_tt,Color,'LineWidth',1.5);
% set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% grid on;
% xlabel ('pump power (W)');
% ylabel ('signal gain (dB)');

 % Graphiques

%--------------------MINE--------------------------------------------
% [XX,YY] = meshgrid([0:dz:Lf-dz],[1:1:jj]);

% figure(1)
% mesh(XX,YY,Pp_evol)
% xlabel ('Length of fiber (m)');
% ylabel ('jj (marker of time evolution)','Rotation',-30);
% zlabel ('Evolution of Pump Output Power (W)');

% figure(2)
% mesh(XX,YY,Ps_evol)
% xlabel ('Length of fiber (m)');
% ylabel ('jj (marker of time evolution)','Rotation',30);
% zlabel ('Evolution of Signal Output Power (W)');

% figure(3)
% mesh(XX,YY,N2_evol)
% figure(4)
% mesh(XX,YY,Ap_evol)
% hold on
% mesh(XX,YY,Am_evol)

% figure(5)

% subplot(2,2,1)
% plot([0:dz:Lf-dz],Pp_evol(jj,:))
% grid on
% xlabel ('Position (m)');
% ylabel ('Pump Power (W)');

% subplot(2,2,2)
% plot([0:dz:Lf-dz],Ps_evol(jj,:))
% grid on
% xlabel ('Position (m)');
% ylabel ('Signal Power (W)');

% subplot(2,2,3)
% plot([0:dz:Lf-dz],Am_evol(jj,:))
% hold on
% plot([0:dz:Lf-dz],Ap_evol(jj,:),'r')
% grid on
% xlabel ('Position (m)');
% ylabel ('ASE Power (W)');
% legend('backward', 'forward')

% subplot(2,2,4)
% plot([0:dz:Lf-dz],N2_evol(jj,:)/N)
% grid on
% xlabel ('Position (m)');
% ylabel ('Upper State Population');

% figure(6)
% plot([0:dz:Lf-dz],Gain_evol(jj,:))
% grid on
% xlabel ('Position (m)');
% ylabel ('Signal Gain (dB)');
%-----------------------END MINE---------------------------------------
%=====================================================================

% Pre- treatment of data (help to save some memory space by clearing Am and
% Ap (occupy lots of memory)
ApL_lbd=Ap(:,ndz);               % (W) Forward ASE (lbd) at z=L
Am0_lbd=Am(:,1);                 % (W) Backward ASE (lbd) at z=0
Ap_z=sum(Ap,1);                  % (W) Total Forward ASE (z)  
Am_z=sum(Am,1);                  % (W) Total Backward ASE (z)  
[ind_lbd]=findech(lambda,lbd_s); % index of lbd_s
ApL_lbd_s=ApL_lbd(ind_lbd,1);    % (W) Forward ASE (lbd_s) at z=L

Filter = logical([zeros(1,jj),ones(1,ntot-jj)]);

PsIN(Filter) = [];
PpIN(Filter) = [];
PsL(Filter) = [];
PpL(Filter) = [];
N2([jj+2:ntot],:)=[];            %to save N2 for all jj

% Help to save memory space 
clear k1 k2 k3 k4 vCa vCae Am Ap 

EDFA_RESTAB;

tmp=strcat('.\res_ai\EDFres',int2str(sim_nb));
save(tmp);

clear tt PpIN PsIN PpL PsL ApL ApLs Am0 Lz N2 N1 Pp Ps Am Ap k1 k2 k3 k4 Gain_s Gain_p Gain_evol Pp_evol Ps_evol N2_evol Ap_s_evol Am_s_evol Ap_evol
clear Am_evol Gain_P_evol

