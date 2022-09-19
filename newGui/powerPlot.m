
%Load data here
global file_selec;



%tmp=strcat('.\res_ai\EDFres',int2str(sim_nb));

load(file_selec);


fs=10; % fontsize
cla (handles.axes1,'reset');
axes(handles.axes1);
if(power_unit)
	Color='b';
	Color2='r';
	[AZ,H1,H2]=plotyy([0:dz:Lf-dz],1000*Ps_evol(jj,:),[0:dz:Lf-dz],1000*Pp_evol(jj,:));
	set(H1,'Color',Color,'LineWidth',1.5);
	set(H2,'Color',Color2,'LineWidth',1.5);
	set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	set(AZ(1),'box','off')
	xlabel(AZ(1),'z(m)','Fontsize',fs);
	ylabel(AZ(1),'Signal Power (mW)','Fontsize',fs);
	set(AZ(2),'box','off')
	ylabel(AZ(2),'Pump Power (mW)','Fontsize',fs);
	grid on;

else
    
	Color='b';
	Color2='r';
	[AZ,H1,H2]=plotyy([0:dz:Lf-dz],10*log10(1000*Ps_evol(jj,:)),[0:dz:Lf-dz],10*log10(1000*Pp_evol(jj,:)));
	set(H1,'Color',Color,'LineWidth',1.5);
	set(H2,'Color',Color2,'LineWidth',1.5);
	set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	set(AZ(1),'box','off')
	xlabel(AZ(1),'z(m)','Fontsize',fs);
	ylabel(AZ(1),'Signal Power (dBm)','Fontsize',fs);
	set(AZ(2),'box','off')
	ylabel(AZ(2),'Pump Power (dBm)','Fontsize',fs);
	grid on;
end



