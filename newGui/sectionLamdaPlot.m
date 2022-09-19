
%Load data here
 tmp = strcat('.\res_ai\EDFres9');
 load(tmp);
 fs=10; % fontsize
cla (handles.axes4,'reset');
axes(handles.axes4);
%====================%%lamda_vs_section========================
        if(choixSection)
		
		    %Clear screen 
			cla 'reset';
            plot(lmbd1,Y1*10^25,'b-','LineWidth',1.5); 
            set(gca,'Fontsize',fs,'FontName','times','XColor','k','YColor','k');
            axis([900 1050 0 3]);
            grid on;
            xlabel ('\lambda [nm]');
            ylabel ('[10^{-25}.m^2]');
		 
		 %====================%%lamda_vs_section2========================
		else
		
		
		Color='b';
	    Color2='r';
		[AZ,H1,H2]=plotyy(lambda*1e9,S_as*10^25,lambda*1e9,S_es*10^25);
		set(H1,'Color',Color,'LineWidth',1.5);
		set(H2,'Color',Color2,'LineWidth',1.5);
		set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
		set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
		set(AZ(1),'box','off')
		xlabel(AZ(1),'\lambda [nm]','Fontsize',fs);
		ylabel(AZ(1),'\sigma_{a} [10^{-25}.m^2]','Fontsize',fs);
		set(AZ(2),'box','off')
		ylabel(AZ(2) , '\sigma_{e} [10^{-25}.m^2]','Fontsize',fs);
		grid on;
		
	    end