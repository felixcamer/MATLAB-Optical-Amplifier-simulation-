
%-----------------ASE POWER--------------------------------------------------------

%Load data here

global file_selec;

%tmp=strcat('.\res_ai\EDFres',int2str(sim_nb));
load(file_selec);
fs=10; % fontsize
cla (handles.axes2,'reset');
axes(handles.axes2);
if(ASE_unit)
	Color='b';
	Color2='r';
	[AZ,H1,H2]=plotyy([0:dz:Lf-dz],1000*Am_evol(jj,:),[0:dz:Lf-dz],1000*Ap_evol(jj,:));
	set(H1,'Color',Color,'LineWidth',1.5);
	set(H2,'Color',Color2,'LineWidth',1.5);
	set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	set(AZ(1),'box','off')
	xlabel(AZ(1),'z(m)','Fontsize',fs);
	ylabel(AZ(1),'Backward ASE Power(mW)','Fontsize',fs);
	set(AZ(2),'box','off')
	ylabel(AZ(2),'Forward ASE Power (mW)','Fontsize',fs);
	grid on;

else
	Color='b';
	Color2='r';
	[AZ,H1,H2]=plotyy([0:dz:Lf-dz],10*log10(1000*Am_evol(jj,:)),[0:dz:Lf-dz],10*log10(1000*Ap_evol(jj,:)));
	set(H1,'Color',Color,'LineWidth',1.5);
	set(H2,'Color',Color2,'LineWidth',1.5);
	set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	set(AZ(1),'box','off')
	xlabel(AZ(1),'z(m)','Fontsize',fs);
	ylabel(AZ(1),'Backward ASE Power(dBm)','Fontsize',fs);
	set(AZ(2),'box','off')
	ylabel(AZ(2),'Forward ASE Power (dBm)','Fontsize',fs);
	grid on;
end