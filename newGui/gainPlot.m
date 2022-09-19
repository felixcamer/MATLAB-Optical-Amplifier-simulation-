%Load data here
 global sim_nb_;
 global power;
 global length_;
 global power_step;
 global length_step;
 global fix_check;
 global start_power_;
 global end_power_;
 
 global end_length_;
 global start_length_;



fs=10; % fontsize
cla (handles.axes3,'reset');
axes(handles.axes3);




%Load data here

tab_gs = [];

tab_gs_hh = [];
len_tab = [];
n = 0;

mod_index = ((end_power_-start_power_)/power_step)+1;

if(fix_check)    %% fix_check = 1   for Ps

	
	for j = 1:sim_nb_

		tmp=strcat('.\res_ai\EDFres',int2str(j));
		load(tmp);
		
		tab_gs = [tab_gs , Gain_evol(end,end)] 
		

		if(mod(j,mod_index)==0)
          n = n+1;
          
          tab_gs_hh(n,:) = tab_gs + 30;
		  len_tab = [len_tab Lf];
		  tab_gs = [];
		 
		 
		end
    end
		
		
	
	 % This script plots force vs displacement stored in the Data struct. 
			%figure, 
			hold on
		    h = legend('show','location','best');
			 set(h,'FontSize',8);
			
			 for u=1:end_length_ -start_length_+ 1 % Data is a struct containing my data.
					
				 name = strcat('Lf= ',int2str(len_tab(u)),'m');
				 plot(start_power_:power_step:end_power_,tab_gs_hh(u,:),'-','LineWidth',1.5,'MarkerSize',1.5,'DisplayName',name);
				
			 end
		    xlabel('PsIN(mw)');
		    ylabel('Gs(dBm)');
	        grid on;
		    
else     %% fix_check = 0   for Pp

       
       
	for j = 1:sim_nb_

		tmp=strcat('.\res_ai\EDFres',int2str(j));
		load(tmp);
		
		tab_gs = [tab_gs , Gain_evol(end,end)]; 
		
		if(mod(j,mod_index)==0)
          n = n+1;
          len_tab = [len_tab Lf];
          tab_gs_hh(n,:) = tab_gs + 30;
		
        
		  tab_gs = [];
		 
		 
		end
    end
	
	 % This script plots force vs displacement stored in the Data struct. 
			%figure, 
			 hold on
			 h = legend('show','location','best');
			set(h,'FontSize',10);
			
			 for u= 1:end_length_ -start_length_+ 1  % Data is a struct containing my data.
					
				 name = strcat('Lf= ',int2str(len_tab(u)),'m');
				 plot(start_power_:power_step:end_power_,tab_gs_hh(u,:),'-','LineWidth',1.5,'MarkerSize',1.5,'DisplayName',name);
				
			 end

			 grid on;
			
		    xlabel('PpIN(mw)');
		    ylabel('Gs(dBm)');
		   
end






