function varargout = EDFA_GUI(varargin)


% EDFA_GUI MATLAB code for EDFA_GUI.fig
%      EDFA_GUI, by itself, creates a new EDFA_GUI or raises the existing
%      singleton*.
%
%      H = EDFA_GUI returns the handle to a new EDFA_GUI or the handle to
%      the existing singleton*.
%
%      EDFA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDFA_GUI.M with the given input arguments.
%
%      EDFA_GUI('Property','Value',...) creates a new EDFA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EDFA_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EDFA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EDFA_GUI

% Last Modified by GUIDE v2.5 23-Apr-2022 07:31:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EDFA_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EDFA_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%==================================GLOBAL VALUES===========================================
             
%==========================================================================================
% --- Executes just before EDFA_GUI is made visible.
function EDFA_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EDFA_GUI (see VARARGIN)

% Choose default command line output for EDFA_GUI
handles.output = hObject;



% handles.laxis = axes('parent',hObject,'units','normalized','position',[0 0 1   1],'visible','off');
% lbls = findobj(hObject,'-regexp','tag','latex_*');
% for i=1: length(lbls)
        % l = lbls(i);
        % % Get current text, position and tag
        % set(l,'units','normalized');
        % s = get(l,'string');
        % p = get(l,'position');
        % t = get(l,'tag');
        % % Remove the UICONTROL
        % delete(l);
        % % Replace it with a TEXT object
        % handles.(t) = text(p(1),p(2),s,'interpreter','latex');
        % disp('once')
% end




% Update handles structure
guidata(hObject, handles);
%set(handles.telecom, 'value', 1);

%--------------afficher la valeur defaut Gs,Gp-------------------
set(handles.gains_display, 'String', "Gsig (dB):" + num2str(0.0));
set(handles.gainp_display, 'String', "Gpump(dB):  " + num2str(0.0));
set(handles.PsdBm, 'String', "Psig(L)(dBm):  " + num2str(0.0));
set(findall(handles.pumpPowermW, '-property', 'enable'), 'enable', 'off');
 
 
set(findall(handles.co_pump, '-property', 'enable'), 'enable', 'off')
set(findall(handles.step_length_edit, '-property', 'enable'), 'enable', 'on')
%set(findall(handles.step_ps_edit, '-property', 'enable'), 'enable', 'on')
%set(findall(handles.step_pp_edit, '-property', 'enable'), 'enable', 'on')


%----test------------------


%---------end---------------



%%% Display the logo.
try
    img = imread('1.png');
    image(handles.fiber_image_axes, img);
    axis(handles.fiber_image_axes, 'off');
    %axis image
	
catch
    axis(handles.fiber_image_axes, 'off');
    warning('img image missing.');
end


%%% Display the logo.
try
    img = imread('ups.png');
    image(handles.axes10, img);
    axis(handles.axes10, 'off');
    %axis image
catch
    axis(handles.axes10, 'off');
    warning('img image missing.');
end

%-----------CLEAR ALL THE SCREENS------------------------------------------
cla (handles.axes1,'reset');
cla (handles.axes2,'reset');
cla (handles.axes3,'reset');
cla (handles.axes4,'reset');

     global sim_nb_;  % this variable is for plotting and save the maximun file number
	 global power;
	 global length_;
	 global power_step;
	 global length_step;
	 global fix_check;
	 global start_power_;
	 global end_power_;
	 global end_length_;
	 global start_length_;
	 global file_selec;

 

% set(handles.axes4,'UserData',get(handles.axes4,'Position'))
% set(handles.axes4,'Position',get(handles.axes4,'UserData'))
%newfig = copyobj(handles.axes4,0); % copy the the figure to an exact new one
%set(handles.axes4, 'Resize', 'on'); % the copy is resizable
%==============thake the freq of pump and signal ==========================

lamdaPp = str2num(get(handles.lamdaPump,'String'));
freg_pump =  round((3*1e8/(lamdaPp*1e-9))/1e12);
%set(handles.freq_pump,freg_pump );
set(handles.freq_pump, 'String', freg_pump);

lamdaPs = str2num(get(handles.lamdaSignal,'String'));
freg_signal =  round((3*1e8/(lamdaPs*1e-9))/1e12);
set(handles.freqSignal, 'String', freg_signal);

%========================END==============================================


% UIWAIT makes EDFA_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = EDFA_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pump_powermW_Callback(hObject, eventdata, handles)
% hObject    handle to pump_powermW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pump_powermW as text
%        str2double(get(hObject,'String')) returns contents of pump_powermW as a double


% --- Executes during object creation, after setting all properties.
function pump_powermW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pump_powermW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function power_indBm_Callback(hObject, eventdata, handles)
% hObject    handle to power_indBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of power_indBm as text
%        str2double(get(hObject,'String')) returns contents of power_indBm as a double


% --- Executes during object creation, after setting all properties.
function power_indBm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to power_indBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_ed_Callback(hObject, eventdata, handles)
% hObject    handle to length_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_ed as text
%        str2double(get(hObject,'String')) returns contents of length_ed as a double


% --- Executes during object creation, after setting all properties.
function length_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pump_wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to pump_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pump_wavelength as text
%        str2double(get(hObject,'String')) returns contents of pump_wavelength as a double


% --- Executes during object creation, after setting all properties.
function pump_wavelength_CreateFcn(hObject, eventdata, ~)
% hObject    handle to pump_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function signal_wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to signal_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of signal_wavelength as text
%        str2double(get(hObject,'String')) returns contents of signal_wavelength as a double


% --- Executes during object creation, after setting all properties.
function signal_wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to signal_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_algo.
function run_algo_Callback(hObject, eventdata, handles)


     %LOCAL VARIABLE INSIDE THE RUN BUTTON
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
     global file_selec; 

      %POWER
      Ps               = str2num(get(handles.powerSignaldBm,'String'));
	  Pp               = str2num(get(handles.pumpPowermW,'String'));
	  start_power      = str2num(get(handles.start_power,'String'));
	  end_power        = str2num(get(handles.end_power,'String'));
	  step_power       = str2num(get(handles.step_power,'String'));
	  
	  power_step = step_power;
	  start_power_ = start_power;
	  end_power_   = end_power;
	  %LAMBDA
	  lamdaPp = str2num(get(handles.lamdaPump,'String'));
	  lamdaPs = str2num(get(handles.lamdaSignal,'String'));
	  
	  
	  
	  %%LENGTHS
	  step_length = str2num(get(handles.step_length_edit,'String'));
	  start_length = str2num(get(handles.start_length,'String'));
	  end_length = str2num(get(handles.end_length,'String'));
	  
	  
	  
	  %VALEUR PAR DEFAULT
	  %step_power = 10;
      length_step = 1;   
      fix_check = 1;
	  %%%%%%%%%%%%%%%%%%%%%%%
	  length_ = end_length;
      
	  end_length_ = end_length;
	  start_length_ = start_length;

	  
	  %check which power is stepped
	  % if(get(handles.step_pp_check,'Value')==0)
	       
		   % disp('One');
		   % fix_check = 0;
		   % power_step  = str2num(get(handles.step_pp_edit,'String'));
		   % power       = Pp;
		   
		
	  % elseif (get(handles.step_ps_check,'Value')==0)
	  
	        % disp('Two');
	       % fix_check = 1;                %the ps is not fixed
		   % power_step  = str2num(get(handles.step_ps_edit,'String'));
		   % power       = Ps;
	  % end
	  
	  
	  %%Toggle buttons (Pp and Ps)
	   if(get(handles.Pp_en,'Value')==1)
	       
		   fix_check = 0;
		   step_power  = str2num(get(handles.step_power,'String'));
		   power       = Pp;
		   
		   set(findall(handles.powerSignaldBm, '-property', 'enable'), 'enable', 'on');
		   set(findall(handles.pumpPowermW, '-property', 'enable'), 'enable', 'off');
		   
	  elseif (get(handles.Ps_en,'Value')==1)
	  
	       fix_check = 1;                %the ps is not fixed
		   step_power  = str2num(get(handles.step_power,'String'));
		   set(findall(handles.pumpPowermW, '-property', 'enable'), 'enable', 'on');
		   set(findall(handles.powerSignaldBm, '-property', 'enable'), 'enable', 'off');
		   power       = Ps;
	  end
	  
	  
	
	%Check the figure type 
	switch get(handles.listFiberFigureType,'Value')


    case 1  %type pump one
		  
		  
	    %-----------CLEAR ALL THE SCREENS----------------------
		 choixDeFibre = (get(handles.TypeFiber,'Value'));
[sim_nb]=  Main_EDFA((choixDeFibre-1),start_power,end_power,Ps,Pp,start_length,end_length,step_length,lamdaPp,lamdaPs, step_power,fix_check);
	sim_nb_ = sim_nb;	 
	file_selec ='.\res_ai\EDFres1' ;
	case 2  %type pump two
	      disp('fig 2 ');
				
	case 3   %type pump tree
	      disp('fig 3 ');
				 
	otherwise %default figure to be figure 1
			
				disp('figure 1');
	end	

	
	%----------WAIT THE PROGRESS BAR-------------	
	waitBar;
	
	%----------PLOT AFTER PROCESSIN--------------
	
	        power_unit = 0;
			ASE_unit = 1;
			gain_unit = 1;
			ASEPlot
			powerPlot
			gainPlot
			N_unit = 1;
			NPlot;
			
			set(handles.gains_display, 'String', "Gsig(dB):   " + num2str(Gain_s));
			set(handles.gainp_display, 'String', "Gpump(dB):  " + num2str(Gain_p));
			set(handles.PsdBm, 'String', "Psig(L)(dBm):"      + num2str(10*log10(1000*PsL(1,jj))));
			set(handles.APLZ, 'String', "Pasefor(L)(dBm)+ :   " + num2str(10*log10(1000*ApL(1,jj))));
			
			F = (2*lbd_s*lbd_s*lbd_s/(h*c*c*dlbd)) * (Ps0(1,jj)* ApLs(1,jj)/(PsL(1,jj)-ApLs(1,jj)));
			set(handles.NF, 'String', "NF(dB):  "      + num2str(10*log10(F)));
			
	%------------END-----------------------------


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in power_length.
function power_length_Callback(hObject, eventdata, handles)
% hObject    handle to power_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of power_length

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% %FOR TELECOM SELECTION
% % --- Executes during object creation, after setting all properties.
% function telecom_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to telecom (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called


% % --- Executes on button press in telecom.
% function telecom_Callback(hObject, eventdata, handles)


% % --- Executes on button press in radix.
% function radix_Callback(hObject, eventdata, handles)


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in relaxation_method.
function relaxation_method_Callback(hObject, eventdata, handles)
% hObject    handle to relaxation_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of relaxation_method



function lamda_Pp_Callback(hObject, eventdata, handles)
% hObject    handle to lamda_Pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lamda_Pp as text
%        str2double(get(hObject,'String')) returns contents of lamda_Pp as a double


% --- Executes during object creation, after setting all properties.
function lamda_Pp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamda_Pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lamda_Ps_Callback(hObject, eventdata, handles)
% hObject    handle to lamda_Ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lamda_Ps as text
%        str2double(get(hObject,'String')) returns contents of lamda_Ps as a double


% --- Executes during object creation, after setting all properties.
function lamda_Ps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamda_Ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in power_vs_time.
function power_vs_time_Callback(hObject, eventdata, handles)
% hObject    handle to power_vs_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of power_vs_time


% --- Executes on button press in nf_vs_ppin.
function nf_vs_ppin_Callback(hObject, eventdata, handles)
% hObject    handle to nf_vs_ppin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nf_vs_ppin


% --- Executes on button press in nf_vs_psin.
function nf_vs_psin_Callback(hObject, eventdata, handles)
% hObject    handle to nf_vs_psin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nf_vs_psin


% --- Executes on button press in gs_vs_ppin.
function gs_vs_ppin_Callback(hObject, eventdata, handles)
% hObject    handle to gs_vs_ppin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gs_vs_ppin


% --- Executes on button press in gs_vs_psin.
function gs_vs_psin_Callback(hObject, eventdata, handles)
% hObject    handle to gs_vs_psin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gs_vs_psin


% --- Executes on button press in ppin_vs_time.
function ppin_vs_time_Callback(hObject, eventdata, handles)
% hObject    handle to ppin_vs_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ppin_vs_time


% --- Executes on button press in lf_vs_n.
function lf_vs_n_Callback(hObject, eventdata, handles)
% hObject    handle to lf_vs_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lf_vs_n


% --- Executes on button press in gain_vs_ps.
function gain_vs_ps_Callback(hObject, eventdata, handles)
% hObject    handle to gain_vs_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gain_vs_ps


% --- Executes on button press in ppin_vs_position.
function ppin_vs_position_Callback(hObject, eventdata, handles)
% hObject    handle to ppin_vs_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ppin_vs_position


% --- Executes on button press in psin_vs_position.
function psin_vs_position_Callback(hObject, eventdata, handles)
% hObject    handle to psin_vs_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of psin_vs_position


% --- Executes on button press in n2_ions_par_m3_vs_posi.
function n2_ions_par_m3_vs_posi_Callback(hObject, eventdata, handles)
% hObject    handle to n2_ions_par_m3_vs_posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of n2_ions_par_m3_vs_posi


% --- Executes on button press in ase_vs_posi.
function ase_vs_posi_Callback(hObject, eventdata, handles)
% hObject    handle to ase_vs_posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ase_vs_posi


% --- Executes on button press in lamda_vs_section.
function lamda_vs_section_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of lamda_vs_section


% --- Executes on button press in lamda_vs_section_2.
function lamda_vs_section_2_Callback(hObject, eventdata, handles)
% hObject    handle to lamda_vs_section_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lamda_vs_section_2




% --- Executes during object creation, after setting all properties.
function lamda_vs_section_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamda_vs_section_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function power_vs_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to power_vs_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ppin_vs_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppin_vs_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gs_vs_ppin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_vs_ppin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ppin_vs_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppin_vs_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function psin_vs_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psin_vs_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function n2_ions_par_m3_vs_posi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2_ions_par_m3_vs_posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ase_vs_posi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ase_vs_posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listFiberFigureType.
function listFiberFigureType_Callback(hObject, eventdata, handles)
% hObject    handle to listFiberFigureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	 cla (handles.axes1,'reset');set(handles.axes1,'Visible','on');
	 cla (handles.axes2,'reset');set(handles.axes2,'Visible','on');
	 cla (handles.axes3,'reset');set(handles.axes3,'Visible','on');
	 cla (handles.axes4,'reset');set(handles.axes4,'Visible','on');
set(findall(handles.co_pump, '-property', 'enable'), 'enable', 'off')
%%% Display 
switch get(handles.listFiberFigureType,'Value')

	  case 1
	  %--------------FIRST CLEAR THE SCREEN-------------------------------------
					 
	        %set(findall(handles.co_pump, '-property', 'enable'), 'enable', 'off')
			try
				Logo = imread('1.png');
				image(handles.fiber_image_axes, Logo);
				axis(handles.fiber_image_axes, 'off');
				%axis image
				
			catch
				axis(handles.fiber_image_axes, 'off');
				warning('Logo image missing.');
			end

	  case 2
	  
	    %--------------FIRST CLEAR THE SCREEN-------------------------------------
					
			
	        %set(findall(handles.co_pump, '-property', 'enable'), 'enable', 'off')
			try
				Logo = imread('2.png');
				image(handles.fiber_image_axes, Logo);
				axis(handles.fiber_image_axes, 'off');
				%axis image
				%set visibility = 1
				
			catch
				axis(handles.fiber_image_axes, 'off');
				warning('Logo image missing.');
			end
			
	 case 3
	   %--------------FIRST CLEAR THE SCREEN-------------------------------------
					
					 
					 
	    set(findall(handles.co_pump, '-property', 'enable'), 'enable', 'on')
		try
			Logo = imread('3.png');
			image(handles.fiber_image_axes, Logo);
			axis(handles.fiber_image_axes, 'off');
			%axis image
		catch
			axis(handles.fiber_image_axes, 'off');
			warning('Logo image missing.');
		end	
			
	otherwise
	
	  %--------------FIRST CLEAR THE SCREEN-------------------------------------
				
		
		try
			Logo = imread('1.png');
			image(handles.fiber_image_axes, Logo);
			axis(handles.fiber_image_axes, 'off');
			%axis image
		catch
			axis(handles.fiber_image_axes, 'off');
			warning('Logo image missing.');
		end		
		
	end		

% Hints: contents = cellstr(get(hObject,'String')) returns listFiberFigureType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listFiberFigureType


% --- Executes during object creation, after setting all properties.
function listFiberFigureType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listFiberFigureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function powerSignaldBm_Callback(hObject, eventdata, handles)
% hObject    handle to powerSignaldBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of powerSignaldBm as text
%        str2double(get(hObject,'String')) returns contents of powerSignaldBm as a double


% --- Executes during object creation, after setting all properties.
function powerSignaldBm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to powerSignaldBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to powerSignaldBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of powerSignaldBm as text
%        str2double(get(hObject,'String')) returns contents of powerSignaldBm as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to powerSignaldBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pumpPowermW_Callback(hObject, eventdata, handles)
% hObject    handle to pumpPowermW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pumpPowermW as text
%        str2double(get(hObject,'String')) returns contents of pumpPowermW as a double


% --- Executes during object creation, after setting all properties.
function pumpPowermW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pumpPowermW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lamdaPump_Callback(hObject, eventdata, handles)
% hObject    handle to lamdaPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lamdaPump as text
%        str2double(get(hObject,'String')) returns contents of lamdaPump as a double


% --- Executes during object creation, after setting all properties.
function lamdaPump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamdaPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in axes3.
function axes3_list_Callback(hObject, eventdata, handles)

switch get(handles.axes3_list,'Value')

	case 1
		
		gainPlot;

	case 2
	       
	   NFPlot;
 	
	otherwise
		
	end	
	
	
	

% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)


%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in axes4_list.
function axes4_list_Callback(hObject, eventdata, handles)

switch get(handles.axes4_list,'Value')

	case 1
		N_unit = 1;
		NPlot;

	case 2
	         powerTimePlote
			
	case 3
	      choixSection = 1;
		  sectionLamdaPlot;
		 
	case 4	  
		   choixSection = 0;
		   sectionLamdaPlot;
			
	otherwise
		
	end		



% --- Executes during object creation, after setting all properties.
function axes4_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dBm_1.
function dBm_1_Callback(hObject, eventdata, handles)

		 power_unit = 0;
		 powerPlot

% --- Executes on button press in mw_1.
function mw_1_Callback(hObject, eventdata, handles)

		 power_unit = 1;
		 powerPlot


% --- Executes on button press in radiobutton20.
function radiobutton20_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton20


% --- Executes on button press in radiobutton19.
function radiobutton19_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton19


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over freqSignal.
function freqSignal_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to freqSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function lamdaSignal_Callback(hObject, eventdata, handles)
% hObject    handle to lamdaSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lamdaSignal as text
%        str2double(get(hObject,'String')) returns contents of lamdaSignal as a double


% --- Executes during object creation, after setting all properties.
function lamdaSignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamdaSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in mw_2.
function mw_2_Callback(hObject, eventdata, handles)
         ASE_unit = 1;
		 ASEPlot


% --- Executes on button press in dBm_2.
function dBm_2_Callback(hObject, eventdata, handles)
         ASE_unit = 0;
		 ASEPlot


% --- Executes on button press in dB_3.
function dB_3_Callback(hObject, eventdata, handles)
         gain_unit = 1;
		 gainPlot


% --- Executes on button press in nat.
function nat_Callback(hObject, eventdata, handles)
         gain_unit = 0;
		 gainPlot


% --- Executes during object creation, after setting all properties.
function axes3_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in axes33_list.
function axes33_list_Callback(hObject, eventdata, handles)
% hObject    handle to axes33_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns axes33_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from axes33_list


% --- Executes during object creation, after setting all properties.
function axes33_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes33_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function gains_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gains_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gainp_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gainp_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in TypeFiber.
function TypeFiber_Callback(hObject, eventdata, handles)
% hObject    handle to TypeFiber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TypeFiber contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeFiber


% --- Executes during object creation, after setting all properties.
function TypeFiber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeFiber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function PsdBm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PsdBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function APLZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to APLZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function NF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes during object creation, after setting all properties.
function Sap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Gp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Sas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Ses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function fiberTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function typeFibreText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typeFibreText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in step_ps_check.
function step_ps_check_Callback(hObject, eventdata, handles)
% hObject    handle to step_ps_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(findall(handles.step_ps_edit, '-property', 'enable'), 'enable', 'off')


% Hint: get(hObject,'Value') returns toggle state of step_ps_check



function step_ps_edit_Callback(hObject, eventdata, handles)
% hObject    handle to step_ps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_ps_edit as text
%        str2double(get(step_ps_edit,'String')) returns contents of step_ps_edit as a double


% --- Executes during object creation, after setting all properties.
function step_ps_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_ps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in step_length_check.
function step_length_check_Callback(hObject, eventdata, handles)
% hObject    handle to step_length_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(findall(handles.step_length_edit, '-property', 'enable'), 'enable', 'off')
% Hint: get(hObject,'Value') returns toggle state of step_length_check
   %disp('checked'); 



function step_length_edit_Callback(hObject, eventdata, handles)
% hObject    handle to step_length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_length_edit as text
%        str2double(get(hObject,'String')) returns contents of step_length_edit as a double


% --- Executes during object creation, after setting all properties.
function step_length_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in step_pp_check.
function step_pp_check_Callback(hObject, eventdata, handles)
% hObject    handle to step_pp_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(findall(handles.step_pp_edit, '-property', 'enable'), 'enable', 'off')

% Hint: get(hObject,'Value') returns toggle state of step_pp_check


% --- Executes on button press in step_pp_edit.
function step_pp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to step_pp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function step_pp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_pp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function end_length_Callback(hObject, eventdata, handles)
% hObject    handle to end_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_length as text
%        str2double(get(hObject,'String')) returns contents of end_length as a double


% --- Executes during object creation, after setting all properties.
function end_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_power_Callback(hObject, eventdata, handles)
% hObject    handle to start_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_power as text
%        str2double(get(hObject,'String')) returns contents of start_power as a double


% --- Executes during object creation, after setting all properties.
function start_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_power_Callback(hObject, eventdata, handles)
% hObject    handle to end_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_power as text
%        str2double(get(hObject,'String')) returns contents of end_power as a double


% --- Executes during object creation, after setting all properties.
function end_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_power_Callback(hObject, eventdata, handles)
% hObject    handle to step_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_power as text
%        str2double(get(hObject,'String')) returns contents of step_power as a double


% --- Executes during object creation, after setting all properties.
function step_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browser_fn.
function browser_fn_Callback(hObject, eventdata, handles)

global file_selec;


[file,path,index] = uigetfile('*.*','Select a single file at a time');


%[file,path] = uigetfile('*.mat');
if isequal(file,0)
   disp('User selected Cancel');
else

  [folder, baseFileNameNoExt, extension] = fileparts(file);
 
end
%-----------------ASE POWER--------------------------------------------------------

tmp=strcat('.\res_ai\',baseFileNameNoExt);
file_selec = tmp;

ASE_unit = 1;
ASEPlot

power_unit = 1;
powerPlot







function start_length_Callback(hObject, eventdata, handles)
% hObject    handle to start_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_length as text
%        str2double(get(hObject,'String')) returns contents of start_length as a double


% --- Executes during object creation, after setting all properties.
function start_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Pp_en_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pp_en (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes during object deletion, before destroying properties.
function Pp_en_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to Pp_en (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Ps_en_CreateFcn(hObject, eventdata, handles)


% --- Executes during object deletion, before destroying properties.
function Ps_en_DeleteFcn(hObject, eventdata, handles)

% --- Executes on button press in Ps_en.
function Ps_en_Callback(hObject, eventdata, handles)
  set(findall(handles.powerSignaldBm, '-property', 'enable'), 'enable', 'off');
  set(findall(handles.pumpPowermW, '-property', 'enable'), 'enable', 'on');
% Hint: get(hObject,'Value') returns toggle state of Ps_en


% --- Executes on button press in Pp_en.
function Pp_en_Callback(hObject, eventdata, handles)

  set(findall(handles.powerSignaldBm, '-property', 'enable'), 'enable', 'on');
  set(findall(handles.pumpPowermW, '-property', 'enable'), 'enable', 'off');


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
