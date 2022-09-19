
close all;
clear all;

% =========================================================================
% Name : Affichage
% Author : Arnaud Fernandez
% Version du : 29/03/2020
%
% Doped fiber : IXblue Rad [qqs m]
% Signal laser : [1530 - 1565 nm]
% Pump laser : around 976 nm 
% 
% Display the numerical results obtained with multiple dopedfiber routines
% dopedfiber_RadIX_multi1
% Compare small signal gain values at different Pump power (-40 dBm) and fixed Signal
% and Fiber length
% Displays : - pump and signal gain through Er3+ doped fiber propagation (time evolution during 
%            simulation too such as steady state regime can be appreciated) 
%            - pump, signal,ASE and N2 spatial evolution during its propagation in Er3+ fiber
%            - ASE power spectrum density  
%            - calculates NF
%
% =========================================================================
close all; clear all;

Sim_start=68;
Sim_stop=117;
nbsim=Sim_stop-Sim_start+1;
inc=0;
for ii=Sim_start:Sim_stop
inc=inc+1;    
    tmp=strcat('test',int2str(ii),'.mat');
    load(tmp)
    
    tabPp_in(1,inc)=Pp_in_W; % (W)
    tabGain_s(1,inc)=10*log10((Ps(1,ndz)*1e3)/(Ps_in_W*1e3)); % (dB)
    tabGain_p(1,inc)=10*log10((Pp(1,ndz)*1e3)/(Pp_in_W*1e3)); % (dB)
      
    
end    

fs=12;
Color='b';
Color2='r';
figure(1)
subplot(211)
plot(tabPp_in,tabGain_s,Color,'LineWidth',1.5);
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('pump power (W)');
ylabel ('signal gain (dB)');
subplot(212)
plot(tabPp_in,tabGain_p,Color2,'LineWidth',1.5);
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('pump power (W)');
ylabel ('pump gain (dB)');

 
% nbstep=size(PsL,2); % (sd) nb of numerical steps
% 
% Gain_s(1,1:nbstep)=10*log10((PsL*1e3)/(Ps_in_W*1e3)); % (dB)
% Gain_p(1,1:nbstep)=10*log10((PpL*1e3)/(Pp_in_W*1e3)); % (dB)
% clear time
% time=linspace(0,nbstep*dt,nbstep); % (s) time simulation vector
% 
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
% subplot(311)
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
% subplot(312)
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
% 
% subplot(313)
% Color='b';
% fs=12;
% plot(Lz,N2*(1e-25),Color,Lz,N*(1e-25),'LineWidth',1.5);
% set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% grid on;
% xlabel ('length (m)');
% ylabel ('N_{2}.10^{25} (m^{-3})');
% 
% figure(4)
% subplot(111)
% Color='b';
% Color2='r';
% fs=12;
% plot(lambda*1e9,10*log10(Ap(:,ndz)*1e3),Color,lambda*1e9,10*log10(Am(:,1)*1e3),Color2,'LineWidth',1.5);
% %plot(lambda*1e9,10*log10(Am(:,1)*1e3),Color2,'LineWidth',1.5);
% set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% grid on;
% xlabel ('lambda (nm)');
% ylabel ('ASE power (dBm/10GHz)');
% legend('forward','backward')
% 
% figure(5)
% subplot(111)
% Color='b';
% Color2='r';
% fs=12;
% plot(time*1e9,PsL*1e3,Color,'LineWidth',1.5);
% %plot(lambda*1e9,10*log10(Am(:,1)*1e3),Color2,'LineWidth',1.5);
% set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
% grid on;
% xlabel ('time (ns)');
% ylabel ('power (mW)');
% %legend('forward','backward')
% 
% fprintf('--------------------------------\n');
% %fprintf('Elapsed rime (s) : %f\n',telapsed);                % (s) elapsed time
% fprintf('Total simulation time (ms) : %f\n',tps_simul*1e3);   
% fprintf('Fiber propagation travel time (ns) : %f\n',tps_fib*1e9);    
% fprintf('\n');
% fprintf('Signal / Pump laser gain (dB) : %f  %f\n',Gain_s(1,ntot),Gain_p(1,ntot));    
% fprintf('\n');
% fprintf('Signal laser INPUT / OUTPUT power (dBm) : %f  %f\n',10*log10(Ps_in_W*1e3),10*log10(PsL(1,ntot)*1e3));    
% fprintf('Pump laser INPUT / OUTPUT power (dBm) : %f  %f\n',10*log10(Pp_in_W*1e3),10*log10(PpL(1,ntot)*1e3));     
% fprintf('Forward / Back. ASE total power (dBm) : %f  %f\n',10*log10(ApL(1,ntot)*1e3),10*log10(Am0(1,ntot)*1e3)); 
% [ind_lbd]=findech(lambda,lbd_s);
% fprintf('\n');
% fprintf('Forward ASE power @ lambda signal (dBm) : %f\n',10*log10((Ap(ind_lbd,ndz)*1e3))); 
% fprintf('--------------------------------\n');
% 
% %save('test00.mat'); % pump 300 mW
% %save('test01.mat'); % pump 100 mW
% 
% 
