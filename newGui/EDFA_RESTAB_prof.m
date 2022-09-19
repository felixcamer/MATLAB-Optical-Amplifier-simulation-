
restab(2,1)=1;
% Selected choice
restab(4,1)=flag_Sae;      % 1 : display of fiber Sae parameters ; 0: don't display
restab(5,1)=flag_RTevo;    % 1 : display power field evolution in the fiber-length ; 0: don't display
restab(6,1)=ASEcalc;       % 1 : ASE included in the model ; 0 : not included
restab(7,1)=choix_fibre;   % 0 : Telecom fiber ; 1 : Rad fiber  
% Fiber physical parameters
restab(9,1)=N;        % (ions.m^(-3))Densité d’ions Er3+
restab(10,1)=Sep;     % (m2)Section efficace d’émission de la pompe @ 980 nm 
restab(11,1)=Gp;      % (sd)Facteur de recouvrement entre le mode de la pompe et l’aire dopée
restab(12,1)=Gs;      % (sd)Facteur de recouvrement entre le mode du signal et l’aire dopée
restab(13,1)=Aeff;    % (m2)Aire de la zone dopée
restab(14,1)=ap;      % (m-1)Pertes linéiques de la pompe
restab(15,1)=as;      % (m-1)Pertes linéiques de la sonde
restab(16,1)=Tau;     % (s)Metastable level lifetime
restab(17,1)=ng;      % (sd)fiber group index
% Frequently modified : signal, pump, fiber
restab(19,1)=Lf;        % (m) doped fiber length
% pump 
restab(21,1)=lbd_p;     % (m) lambda pump
restab(22,1)=Sap;       % (m2) 
restab(23,1)=Pp_in_dBm; % (dBm) input pump power 
restab(24,1)=Pp_in_W;   % (W) input pump power 
% Signal
restab(26,1)=lbd_s;     % (m) input signal wavelength 
restab(27,1)=Ses;       % (m2)
restab(28,1)=Sas;       % (m2)
restab(29,1)=Ps_in_dBm; % (dBm) input signal power
restab(30,1)=Ps_in_W;   % (W) input signal power
% ESA
restab(32,1)=dlbd;           % step spectral de l'ESA delta lambda (m)
restab(33,1)=dnu;            % step spectral de l'ESA (Hz)
% Numerical sampling
restab(35,1)=dz;          % (m) spatial step index 
restab(36,1)=ndz;         % (sd) nb of spatial steps
restab(37,1)=dt;          % (s) temporal step
restab(38,1)=nlbd;        % (sd) nb of spectral slices  
restab(39,1)=tps_fib;     % (s) Maximum fiber propagation time
restab(40,1)=tps_simul;   % (s) Maximum simulation time 
restab(41,1)=ntot;        % (sd) nb of temporal steps
restab(42,1)=Trise;       % (s) rise time for a slow turn-on of input signal and pump laser 
% Tolerance 
restab(44,1)=tol;         % (W/ns) tolerance
restab(45,1)=counter;     % (nd) nb of required consecutive occurence of below tolerance mtelapsed     
restab(46,1)=tth;         % (s) threshold time for evaluation
restab(47,1)=delT;        % (s) delta t
% Results
restab(49,1)=PpL(1,jj);   % (W) outgoing pump 
restab(50,1)=Ps0(1,jj);   % (W) ingoing signal
restab(51,1)=PsL(1,jj);   % (W) outgoing signal
restab(52,1)=ApL(1,jj);   % (W) outgoing forward ASE @ z=L
restab(53,1)=ApLs(1,jj);  % (W) outgoing forward ASE @ z=L and lambda signal
restab(54,1)=Am0(1,jj);   % (W) outgoing backward ASE @ z=0
restab(55,1)=Gain_s;      % signal gain (dB)
restab(56,1)=Gain_p;      % pump gain (dB)
restab(57,1)=telapsed;    % elapsed time (s)



