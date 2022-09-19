%Load data here
 global sim_nb;
 tmp=strcat('.\res_ai\EDFres',int2str(sim_nb));
 load(tmp);
 fs=10; % fontsize
cla (handles.axes4,'reset');
axes(handles.axes4);
if(N_unit)
	Color='b';
	Color2='r';
	[AZ,H1,H2]=plotyy([0:dz:Lf-dz],(N-N2(1,:))*(1e-25),[0:dz:Lf-dz],N2(1,:)*(1e-25));
	set(H1,'Color',Color,'LineWidth',1.5);
	set(H2,'Color',Color2,'LineWidth',1.5);
	set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	set(AZ(1),'box','off')
	xlabel(AZ(1),'z(m)','Fontsize',fs);
	ylabel(AZ(1),'N1 ions/m^3','Fontsize',fs);
	set(AZ(2),'box','off')
	ylabel(AZ(2) , 'N2 ions/m^3','Fontsize',fs);
	grid on;

% else
    
	% Color='b';
	% Color2='r';
	% [AZ,H1,H2]=plotyy([0:dz:Lf-dz],10.^((Gain_evol(jj,:))/20),[0:dz:Lf-dz],10.^((Gain_P_evol(jj,:))/20));
	% set(H1,'Color',Color,'LineWidth',1.5);
	% set(H2,'Color',Color2,'LineWidth',1.5);
	% set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
	% set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
	% set(AZ(1),'box','off')
	% xlabel(AZ(1),'z(m)','Fontsize',fs);
	% ylabel(AZ(1),'Signal Gain (nat)','Fontsize',fs);
	% set(AZ(2),'box','off')
	% ylabel(AZ(2) , 'Pump Gain (nat)','Fontsize',fs);
	% grid on;

end



