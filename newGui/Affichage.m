
close all;
clear all;

% =========================================================================
% Name : Affichage
% Author : Arnaud Fernandez
% Version du : 03/06/2020
%
% Doped fiber : IXblue Rad ou Telecom [qqs m]
% Signal laser : [1530 - 1565 nm]
% Pump laser : around 976 nm 
% 
% Display the numerical results obtained with dopedfiber routines 
% Displays : - pump and signal gain through Er3+ doped fiber propagation (time evolution during 
%            simulation too such as steady state regime can be appreciated) 
%            - pump, signal,ASE and N2 spatial evolution during its propagation in Er3+ fiber
%            - ASE power spectrum density  
%            - calculates NF
%
% =========================================================================


%tmp = strcat('.\res\Resultats de simu_prof_drive\EDFres4.mat');
%tmp = strcat('.\res_prof\EDFres7');
tmp = strcat('.\res_ai\EDFres9');
load(tmp)

Color='b';
Color2='r';
fs=12;
figure(1)

%Tracé de la Puissance signal en fonction de temps

subplot(2,2,1)
plot(time(1,1:jj)*1e9,PsIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
hold on
plot(time(1,1:jj)*1e9,PsL(1,1:jj)*1e3,Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('time (ns)');
ylabel ('power (mW)');
legend('input signal power','output signal power')

subplot(2,2,3)
yyaxis left
plot(time(1,1:jj)*1e9,PsIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
ylabel ('input signal power (mW)');
yyaxis right
plot(time(1,1:jj)*1e9,PsL(1,1:jj)*1e3,Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
set(gca,'Fontsize',fs,'FontName','times','XColor','k');
grid on;
xlabel ('time (ns)');
ylabel ('output signal power (mW)');
legend('input signal power','output signal power')

%Tracé de la Puissance pompe en fonction de temps

subplot(2,2,2)
plot(time(1,1:jj)*1e9,PpIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
hold on
plot(time(1,1:jj)*1e9,PpL(1,1:jj)*1e3,Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('time (ns)');
ylabel ('power (mW)');
legend('input pump power','output pump power')

subplot(2,2,4)
yyaxis left
plot(time(1,1:jj)*1e9,PpIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
ylabel ('input pump power (mW)');
yyaxis right
plot(time(1,1:jj)*1e9,PpL(1,1:jj)*1e3,Color2,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
set(gca,'Fontsize',fs,'FontName','times','XColor','k');
grid on;
xlabel ('time (ns)');
ylabel ('output pump power (mW)');
legend('input pump power','output pump power')


%Tracé du GAIN

Gain_s_tt(1,:)=10*log10((PsL*1e3)/(Ps_in_W*1e3));
Gain_p_tt(1,:)=10*log10((PpL*1e3)/(Pp_in_W*1e3));
ntot=jj;

figure(2)
subplot(1,1,1)
[AZ,H1,H2]=plotyy(time(1,1:jj)*1e9,Gain_s_tt,time(1,1:jj)*1e9,Gain_p_tt);
set(H1,'Color',Color,'LineWidth',1.5);
set(H2,'Color',Color2,'LineWidth',1.5);
set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
set(AZ(1),'box','off')
xlabel(AZ(1),'time (ns)','Fontsize',fs);
ylabel(AZ(1),'Signal gain (dB)','Fontsize',fs);
set(AZ(2),'box','off')
ylabel(AZ(2),'Pump gain (dB)','Fontsize',fs);
grid on;

%Tracé du Pase en fonction de la longeur de la fibre

figure(3)
subplot(3,1,1)
plot(Lz,10*log10(Ap_z*1e3),Color,Lz,10*log10(Am_z*1e3),Color2,'LineWidth',1.5);
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('length (m)');
ylabel ('ASE power (dBm)');
legend('forward','backward')

%Tracé de les Puissance signal et pompe en fonction de la longueur de la fibre

subplot(3,1,2) 
Color='b';
Color2='r';
[AZ,H1,H2]=plotyy(Lz,10*log10(Ps*1e3),Lz,10*log10(Pp*1e3));
set(H1,'Color',Color,'LineWidth',1.5);
set(H2,'Color',Color2,'LineWidth',1.5);
set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
set(AZ(1),'box','off')
xlabel(AZ(1),'length (m)','Fontsize',fs);
ylabel(AZ(1),'Signal power (dBm)','Fontsize',fs);
set(AZ(2),'box','off')
ylabel(AZ(2),'Pump power(dBm)','Fontsize',fs);
grid on;


%Tracé du N2 en fonction de la longeur de la fibre

subplot(3,1,3) 
plot(Lz,N2(1,:)*(1e-25),Color,Lz,N*(1e-25),Color2,'LineWidth',1.5);
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('length (m)');
ylabel ('N_{2}.10^{25} (m^{-3})');

%Tracé la densite scectrale de Pase en fonction de la longueur d'onde

figure(4)
subplot(1,1,1)
plot(lambda*1e9,10*log10(ApL_lbd*1e3),Color,lambda*1e9,10*log10(Am0_lbd*1e3),Color2,'LineWidth',1.5);
hold on;
plot(lbd_s*1e9,10*log10(PsL(1,jj)*1e3),'*');
set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
grid on;
xlabel ('lambda (nm)');
ylabel ('ASE power (dBm/10GHz)');
legend('forward','backward')

%Affichafe des parametres

fprintf('--------------------------------\n');
%fprintf('Elapsed rime (s) : %f\n',telapsed);                % (s) elapsed time
fprintf('Total simulation time (s) : %f\n',telapsed);   
fprintf('Fiber propagation travel time (ns) : %f\n',tps_fib*1e9);    
fprintf('\n');
fprintf('Signal / Pump laser gain (dB) : %f  %f\n',Gain_s,Gain_p);    
fprintf('\n');
fprintf('Signal laser INPUT / OUTPUT power (dBm) : %f  %f\n',10*log10(Ps_in_W*1e3),10*log10(PsL(1,ntot)*1e3));    
fprintf('Pump laser INPUT / OUTPUT power (dBm) : %f  %f\n',10*log10(Pp_in_W*1e3),10*log10(PpL(1,ntot)*1e3));     
fprintf('Forward / Back. ASE total power (dBm) : %f  %f\n',10*log10(ApL(1,ntot)*1e3),10*log10(Am0(1,ntot)*1e3)); 
fprintf('\n');
fprintf('Forward ASE power @ lambda signal (dBm) : %f\n',10*log10((ApL_lbd_s*1e3))); 
fprintf('--------------------------------\n');


%Tracé l'augmentation de la Puissance pompe en fonction de temps

figure(5)
yyaxis left
plot(time(1,1:jj)*1e9,PpIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','--')
ylabel ('input pump power (mW)');
yyaxis right
plot(time(1,1:jj)*1e9,PsL(1,1:jj)*1e3,Color2,time(1,1:jj)*1e9,PsIN(1,1:jj)*1e3,Color,'markerfacecolor',[1 0 0],'LineWidth',2.5,'LineStyle','-')
set(gca,'Fontsize',fs,'FontName','times','XColor','k');
grid on;
xlabel ('time (ns)');
ylabel ('output pump power (mW)');
legend('input pump power','output signal power', 'input signal power')


figure(6)
[X,Y] = meshgrid(Lz,[1:jj]);
mesh(X,Y,N2(2:jj+1,:)*(1e-25));