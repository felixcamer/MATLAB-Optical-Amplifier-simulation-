% close all;
% clear all;

% =========================================================================
% Name : Main_EDFA
% Author : Barkoum Btera Felix
% Version du : 03/04/2022
%
% Doped fiber : IXblue Rad or TelecomIX
% Signal laser : [1530 - 1565 nm]
% Pump laser : around 976 nm 
% 
% Calculates : forward and backward ASE ; cw signal laser amplification ; NF 
% 
% Physical properties integrated by the model :
% cw pump depletion along the doped fiber ; Er3+ ions density dependance over length ;
% Er3+ transition in glass is homogeneously broadened.
% 
% Call : Propag_EDFA that computes EDF propagation equation and rate equation
% eqs (6.3) to (6.7) p 157 ; EDFA [Becker]
%
% =========================================================================


function  [sim_nb ]=  Main_EDFA(choix_fibre,start_power,end_power,Ps_maxmw,Pp_maxmw,start_length,end_length,step_length,lamdaPp,lamdaPs, step_power,fix_check)

flag_Sae=0;         % 1 : display of fiber Sae parameters ; 0: don't display
flag_RTevo=0;       % 1 : display power field evolution in the fiber-length ; 0: don't display
ASEcalc=1;          % 1 : ASE included in the model ; 0 : not included
%choix_fibre=1;      % 0 : Telecom fiber ; 1 : Rad fiber  


switch choix_fibre
    
    case 0  
        
        % Load parameters  iXBlueTelecom ----------------------------------------------
	     root='./TelecomIX/data1.txt';
         A  =  load(root);
         lmbd1 = A(:,1);
         Y1 = A(:,2);
         if flag_Sae==1
            figure(1)
            subplot(111)
            fs=20; % fontsize
            plot(lmbd1,Y1*10^25,'b-','LineWidth',1.5); 
            set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
            axis([900 1050 0 2.6]);
            legend('\sigma_{ap}');
            grid on;
            xlabel ('\lambda [nm]');
            ylabel ('[10^{-25}.m^2]');
            
            %title('IX fiber : Rad hard');
         end
		 
          Y1 = Y1';    %THIS IS FOR TELECOM FIBER  IF NOT AN ERROR OCCUR
		  
		    root='./TelecomIX/iXBlueTelpart2.txt';
			[lambda,S_as,S_es]=my_data4(root); % [nm,m2,m2] lambda,Sas,Ses
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
		 
		 
        % [lambda,S_as,S_es]=my_data4(root); % [nm,m2,m2] lambda,Sas,Ses
        % lambda=lambda*1e-9; % (m)

         % if flag_Sae==1
            % figure(2)
            % subplot(111)
            % fs=20; % fontsize
            % plot(lambda*1e9,S_as*10^25,'b-',lambda*1e9,S_es*10^25,'r-','LineWidth',1.5); 
            % set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
            % axis([1450 1650 0  inf]);
            % legend('\sigma_{a}','\sigma_{e}');
            % grid on;
            % xlabel ('\lambda [nm]');
            % ylabel ('[10^{-25}.m^2]');
            % title('IX fiber : Telecom');
          % end

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
        %S_as = S_as/1.6;
        %S_es = S_es/1.6;
        
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
		
		%-------------------END RAD FIBER---------------------------
		
		
		case 2   

        % Load parameters  Light-125/8 ----------------------------------------------

  	         root='./nLIGTH/nLIGTH_2.txt';
			 A  =  load(root);
			 lmbd1 = A(:,1);
			 Y1 = A(:,2);
		  
	     if flag_Sae==1
           
			figure(1)
			subplot(111)
			fs=20; % fontsize
			plot(lmbd1,Y1*10^25,'b-','LineWidth',1.5); 
			
			set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
			axis([900 1050 0 2.9]);
			legend('\sigma_{ap}');
			grid on;
			xlabel ('\lambda [nm]');
			ylabel ('[10^{-25}.m^2]');
	        title('nLIGTH_8_125_1');
		  end 
		  
		   Y1 = Y1'; 
		   root='./nLIGTH/nLIGHT_1.txt';
		   [lambda,S_as,S_es]=my_data4(root); % [nm,m2,m2] lambda,Sas,Ses
		   lambda=lambda*1e-9; % (nm)
		   
		   if flag_Sae==1
		
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
            title('nLIGTH_8_125_1');
		   
		  end
		  lambda=lambda*1e-9;
		  
        %%CONSTANTS HERE
			N=3.5e25;       % (ions.m^(-3))	Densit?d’ions Er3+
			Sep=0;          % (m2)	Section efficace d’émission de la pompe @ 980 nm 
			Gp=0.76346;        % (sd) Facteur de recouvrement entre le mode de la pompe et l’aire dopée
			Gs=0.75621;        % (sd) Facteur de recouvrement entre le mode du signal et l’aire dopée
			Aeff=6.4e-23; % (m2)	Aire de la zone dopée
			ap=8.4962e-3;   % (m-1)	Pertes linéiques de la pompe  1.3816e-3  8.5  8.1
			as=8.0756e-3 ;   % (m-1)	Pertes linéiques de la sonde
			Tau=9e-3;      % (s) Metastable level lifetime
			ng=1.5;         % (sd) fiber group index
		

		case 3

         % Load parameters  Light-125/4 ----------------------------------------------

  	         root='./nLIGTH/nLIGTH_2.txt';
			 A  =  load(root);
			 lmbd1 = A(:,1);
			 Y1 = A(:,2);
		 if flag_Sae==1
           
			figure(1)
			subplot(111)
			fs=20; % fontsize
			plot(lmbd1,Y1*10^25,'b-','LineWidth',1.5); 
			
			set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
			axis([900 1050 0 2.9]);
			legend('\sigma_{ap}');
			grid on;
			xlabel ('\lambda [nm]');
			ylabel ('[10^{-25}.m^2]');
	        title('nLIGTH_8_125_1');
		  end 
		  
		   Y1 = Y1'; 
		   root='./nLIGTH/nLIGHT_1.txt';
		   [lambda,S_as,S_es]=my_data4(root); % [nm,m2,m2] lambda,Sas,Ses
		   lambda=lambda*1e-9; % (nm)
		   
		   if flag_Sae==1
		
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
            title('nLIGTH_8_125_1');
		   
		  end
		  lambda=lambda*1e-9;
		  
        %%CONSTANTS HERE
			N=8e25;       % (ions.m^(-3))	Densit?d’ions Er3+
			Sep=0;          % (m2)	Section efficace d’émission de la pompe @ 980 nm 
			Gp=0.45523;        % (sd) Facteur de recouvrement entre le mode de la pompe et l’aire dopée
			Gs=0.44333;        % (sd) Facteur de recouvrement entre le mode du signal et l’aire dopée
			Aeff=2.95e-23; % (m2)	Aire de la zone dopée
			ap=15.412e-3;   % (m-1)	Pertes linéiques de la pompe
			as=11.516e-3;   % (m-1)	Pertes linéiques de la sonde
			Tau=9e-3;      % (s) Metastable level lifetime
			ng=1.5;         % (sd) fiber group index
		
		
		
		case 4

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
        %S_as = S_as/1.6;
        %S_es = S_es/1.6;
        
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

			N=1.4e25;       % (ions.m^(-3))	Densit?d’ions Er3+
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

%==========================================================================
% Initial conditions

lbd_p=lamdaPp*1e-9;                          % (m) lambda pump
nu_p=c/lbd_p;                            % (Hz) pump carrier frequency 
[Sap] = findSap(lbd_p,lmbd1*1e-9,Y1);    % (m2) 
lbd_s=lamdaPs*1e-9;                           % (m) input signal wavelength 
nu_s=c/lbd_s;                            % (Hz) signal carrier frequency 
dlbd=abs(lambda(2)-lambda(1));           % delta lambda (m)
dnu=(c/(lbd_s^2))*dlbd;                  % (Hz) delta nu
nu=c./lambda;                            % (Hz) freq. vector


% intermediate calculations for propagation equations

[Sas,Ses]=findSasSes(lbd_s,lambda,S_as,S_es);     % (m2)
Cas=(Tau*Sas*Gs)/(Aeff*h*nu_s);
Cap=(Tau*Sap*Gp)/(Aeff*h*nu_p); 
C=(Tau*Gs)/(Aeff*h);
vCa=(lambda'/c).*(S_as)';
vCae=(lambda'/c).*(S_as+S_es)';
Ces=(Tau*Ses*Gs)/(Aeff*h*nu_s);
Cep=(Tau*Sep*Gp)/(Aeff*h*nu_p);         

% --------------------------------------------------------------
% RES : xx to xx
%Modifiable
% Ps_in_dBm=-Psin_dBm;                    % (dBm) input signal power -40
% Pp_in_W= Pp_in_mW*0.001;                % (W) input pump power  140e-3
% Ps_in_W=(10.^(Ps_in_dBm/10))*(1e-3);    % (W) input signal power
% Pp_in_dBm=10*log10(Pp_in_W*1e3);        % (dBm) pump power
% dz=0.2;                               % (m) spatial step index 
% Lfstart = Length;                             % fiberlength (m)
% Lfstop =  Length;                             % fiberlength (m)
% sinb = 9;
% nbloop = Lfstop - Lfstart + 1;          % number of loops
% PpINstart = Pp_in_W; % en W
% PpINstop = Pp_in_W + 0.01; % en W
% PpINstep = (PpINstop*1000 - PpINstart*1000)/(10) + 1; 
% %initialisation de restab
% restab = zeros(57,PpINstep);
% Ltab=linspace(Lfstart,Lfstop,nbloop);
% PpINtab=linspace(PpINstart,PpINstop,PpINstep);
% for j = 1 : 1 : nbloop
   % Lf=Ltab(1,j);
   % fprintf('Simu. : %f\n',sim_nb);
   % Propag_EDFA;
   % sim_nb = sim_nb + 1;
   % %pause;
   % % Need to declare vCa and vCae  that were cleared in Propag_EDFA before saving 
   % vCa=(lambda'/c).*(S_as)';
   % vCae=(lambda'/c).*(S_as+S_es)';
% end
% --------------------------------------------------------------
% Ps_in_dBm=-Psin_dBm;                    % (dBm) input signal power -40
% Ps_in_W=(10.^(Ps_in_dBm/10))*(1e-3);    % (W) input signal power


if(fix_check)  %%select the Ps

	dz=0.2;% (m) spatial step index 
	%lenght here
	sim_nb = 0;
	nb_lenght_step = step_length;

	%Pump Input signal 
	Pp_in_W   =  Pp_maxmw*1e-3;% (W) input pump power  140e-3  Pp_in_mW*0.001; 
	Pp_in_dBm =  10*log10(Pp_in_W*1e3); % (dBm) pump power
	PpINstart =  Pp_in_W; % en W
	PpINstop  =  Pp_in_W + 0.01; %en W
	PpINstep  = (PpINstop*1000 - PpINstart*1000)/(10) + 1; 
	PpINtab   =  linspace(PpINstart,PpINstop,PpINstep);
	restab    = zeros(57,PpINstep);
	for Lf = start_length:nb_lenght_step: end_length
	   
      
	   for i = start_power:step_power:end_power
		Ps_in_W = i*1e-3;
		Ps_in_dBm=10*log10(i); 
		sim_nb = sim_nb + 1
		%Modifiable
	    Propag_EDFA;
	    vCa=(lambda'/c).*(S_as)';
	    vCae=(lambda'/c).*(S_as+S_es)';
	   end
	 end
	 
else
	 
	 
	dz=0.2;% (m) spatial step index 
	sim_nb = 0;
	nb_lenght_step = step_length;

	%Pump Input signal
	Ps_in_W = Ps_maxmw*1e-3;
	Ps_in_dBm=10*log10(Ps_in_W); 
	
	
	 for Lf = start_length:nb_lenght_step: end_length
	   
      
	   for i = start_power:step_power:end_power
	   
	   Pp_in_W   = i*1e-3;% (W) input pump power  140e-3  Pp_in_mW*0.001; 
	   Pp_in_dBm =  10*log10(Pp_in_W*1e3); % (dBm) pump power
	   PpINstart =  Pp_in_W; % en W
	   PpINstop  =  Pp_in_W + 0.01; %en W
	   PpINstep  = (PpINstop*1000 - PpINstart*1000)/(10) + 1; 
	   PpINtab   =  linspace(PpINstart,PpINstop,PpINstep);
	   restab    = zeros(57,round(PpINstep));
		
		sim_nb = sim_nb + 1;
		
		%Modifiable
	    Propag_EDFA;
	    vCa=(lambda'/c).*(S_as)';
	    vCae=(lambda'/c).*(S_as+S_es)';
	   end
	 end
	
 end
