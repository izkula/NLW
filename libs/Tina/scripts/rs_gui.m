function varargout = rs_gui(varargin)
% RS_GUI MATLAB code for rs_gui.fig
%      RS_GUI, by itself, creates a new RS_GUI or raises the existing
%      singleton*.
%
%      H = RS_GUI returns the handle to a new RS_GUI or the handle to
%      the existing singleton*.
%
%      RS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RS_GUI.M with the given input arguments.
%
%      RS_GUI('Property','Value',...) creates a new RS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rs_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rs_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rs_gui

% Last Modified by GUIDE v2.5 15-Dec-2016 14:17:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rs_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @rs_gui_OutputFcn, ...
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


% --- Executes just before rs_gui is made visible.
function rs_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rs_gui (see VARARGIN)

% Choose default command line output for rs_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rs_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rs_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function enter_pathname_Callback(hObject, eventdata, handles)
% hObject    handle to enter_pathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_pathname as text
%        str2double(get(hObject,'String')) returns contents of enter_pathname as a double
pathname=get(hObject,'String');
display(pathname);
handles.pathname = pathname;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_pathname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_pathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.pathname = 'F:/Users/Tina/Dropbox/MATLAB/DEISSEROTH/rs_scope/data/';
guidata(hObject,handles)


function enter_filename_Callback(hObject, eventdata, handles)
% hObject    handle to enter_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_filename as text
%        str2double(get(hObject,'String')) returns contents of enter_filename as a double
filename=get(hObject,'String');
handles.filename = filename;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.settings.filename = '';
guidata(hObject,handles)


function enter_filedate_Callback(hObject, eventdata, handles)
% hObject    handle to enter_filedate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_filedate as text
%        str2double(get(hObject,'String')) returns contents of enter_filedate as a double
filedate=get(hObject,'String');
handles.filedate = filedate;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_filedate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_filedate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.filedate = '';
guidata(hObject,handles)


function enter_flims_Callback(hObject, eventdata, handles)
% hObject    handle to enter_flims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_flims as text
%        str2double(get(hObject,'String')) returns contents of enter_flims as a double
flims=get(hObject,'String');
if ~isempty(flims)
    flims_split=strsplit(flims,':');
    handles.flims = [str2double(flims_split{1}):str2double(flims_split{2})];
else
    handles.flims=[];
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_flims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_flims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.flims = [];
guidata(hObject,handles)


function enter_nPCs_Callback(hObject, eventdata, handles)
% hObject    handle to enter_nPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_nPCs as text
%        str2double(get(hObject,'String')) returns contents of enter_nPCs as a double
nPCs=str2double(get(hObject,'String'));
handles.nPCs=nPCs;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_nPCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_nPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.nPCs=200;
guidata(hObject,handles)


function enter_dsamp_Callback(hObject, eventdata, handles)
% hObject    handle to enter_dsamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_dsamp as text
%        str2double(get(hObject,'String')) returns contents of enter_dsamp as a double
dsamp=get(hObject,'String');
if ~isempty(dsamp)
    dsamp_split=strsplit(dsamp,':');
    handles.dsamp = [str2double(dsamp_split{1}):str2double(dsamp_split{2})];
else
    handles.dsamp=[];
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_dsamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_dsamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dsamp=[];
guidata(hObject,handles)


function enter_badframes_Callback(hObject, eventdata, handles)
% hObject    handle to enter_badframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_badframes as text
%        str2double(get(hObject,'String')) returns contents of enter_badframes as a double
badframes=get(hObject,'String');
if ~isempty(badframes)
    badframes_split=strsplit(badframes,':');
    badframes.dsamp = [str2double(badframes_split{1}):str2double(badframes_split{2})];
else
    badframes.dsamp = [];
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_badframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_badframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.badframes=[];
guidata(hObject,handles)


function enter_mu_Callback(hObject, eventdata, handles)
% hObject    handle to enter_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_mu as text
%        str2double(get(hObject,'String')) returns contents of enter_mu as a double
mu=str2double(get(hObject,'String'));
handles.mu=mu;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.mu=0.1;
guidata(hObject,handles)

% --- Executes on button press in run_pcaica.
function run_pcaica_Callback(hObject, eventdata, handles)
% hObject    handle to run_pcaica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_pcaica(handles.pathname,handles.filedate,handles.filename,handles.flims,...
    handles.nPCs,handles.dsamp,handles.badframes,handles.mu);


function enter_smwidth_Callback(hObject, eventdata, handles)
% hObject    handle to enter_smwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_smwidth as text
%        str2double(get(hObject,'String')) returns contents of enter_smwidth as a double
smwidth=str2double(get(hObject,'String'));
handles.smwidth=smwidth;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_smwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_smwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.smwidth=0.2;
guidata(hObject,handles)


function enter_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to enter_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_thresh as text
%        str2double(get(hObject,'String')) returns contents of enter_thresh as a double
thresh=str2double(get(hObject,'String'));
handles.thresh=thresh;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.thresh=4;
guidata(hObject,handles)


function enter_arealims_Callback(hObject, eventdata, handles)
% hObject    handle to enter_arealims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_arealims as text
%        str2double(get(hObject,'String')) returns contents of enter_arealims as a double
arealims=str2double(get(hObject,'String'));
handles.arealims=arealims;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_arealims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_arealims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.arealims=20;
guidata(hObject,handles)


function enter_plotting_Callback(hObject, eventdata, handles)
% hObject    handle to enter_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_plotting as text
%        str2double(get(hObject,'String')) returns contents of enter_plotting as a double
plotting=str2double(get(hObject,'String'));
handles.plotting=plotting;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_plotting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.plotting=1;
guidata(hObject,handles)


function enter_frametime_Callback(hObject, eventdata, handles)
% hObject    handle to enter_frametime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_frametime as text
%        str2double(get(hObject,'String')) returns contents of enter_frametime as a double
frametime=str2double(get(hObject,'String'));
handles.frametime=frametime;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_frametime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_frametime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.frametime=0.132423184;
guidata(hObject,handles)


function enter_deconvtauoff_Callback(hObject, eventdata, handles)
% hObject    handle to enter_deconvtauoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_deconvtauoff as text
%        str2double(get(hObject,'String')) returns contents of enter_deconvtauoff as a double
deconvtauoff=str2double(get(hObject,'String'));
handles.deconvtauoff=deconvtauoff;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_deconvtauoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_deconvtauoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.deconvtauoff=0.1;
guidata(hObject,handles)


function enter_slidingwind_Callback(hObject, eventdata, handles)
% hObject    handle to enter_slidingwind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_slidingwind as text
%        str2double(get(hObject,'String')) returns contents of enter_slidingwind as a double
slidingwind=str2double(get(hObject,'String'));
handles.slidingwind=slidingwind;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_slidingwind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_slidingwind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.slidingwind=10;
guidata(hObject,handles)


function ener_percentfilter_Callback(hObject, eventdata, handles)
% hObject    handle to ener_percentfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ener_percentfilter as text
%        str2double(get(hObject,'String')) returns contents of ener_percentfilter as a double
percentfilter=str2double(get(hObject,'String'));
handles.percentfilter=percentfilter;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function ener_percentfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ener_percentfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.percentfilter=8;
guidata(hObject,handles)


function enter_lowpassfilt_Callback(hObject, eventdata, handles)
% hObject    handle to enter_lowpassfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_lowpassfilt as text
%        str2double(get(hObject,'String')) returns contents of enter_lowpassfilt as a double
lowpassfilt=str2double(get(hObject,'String'));
handles.lowpassfilt=lowpassfilt;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function enter_lowpassfilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_lowpassfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.lowpassfilt=3;
guidata(hObject,handles)


% --- Executes on button press in run_applycellmasks.
function run_applycellmasks_Callback(hObject, eventdata, handles)
% hObject    handle to run_applycellmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_applycellmasks(handles.pathname,handles.filedate,handles.filename,handles.frametime,...
    handles.deconvtauoff,handles.slidingwind,handles.percentfilter,...
    handles.lowpassfilt);


% --- Executes on button press in run_plotdfof.
function run_plotdfof_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotdfof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_plotdfof(handles.frametime,handles.pathname,handles.filedate,handles.filename);


% --- Executes on button press in run_deletecellmasks.
function run_deletecellmasks_Callback(hObject, eventdata, handles)
% hObject    handle to run_deletecellmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_deletecellmasks(handles.pathname,handles.filedate,handles.filename);


% --- Executes on button press in run_findcellmasks.
function run_findcellmasks_Callback(hObject, eventdata, handles)
% hObject    handle to run_findcellmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_findcellmasks(handles.pathname,handles.filedate,handles.filename,handles.smwidth,...
    handles.thresh,handles.arealims,handles.plotting);


% --- Executes on button press in run_plotcellmasks.
function run_plotcellmasks_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotcellmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_plotcellmasks(handles.pathname,handles.filedate,handles.filename);


% --- Executes on button press in run_loadbinfile.
function run_loadbinfile_Callback(hObject, eventdata, handles)
% hObject    handle to run_loadbinfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_loadbinfile(handles.pathname,handles.filedate,handles.filename,handles.frametime)


% --- Executes on button press in run_plotleverpress.
function run_plotleverpress_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotleverpress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in run_plotwaterlicking.
function run_plotwaterlicking_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotwaterlicking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in run_plotshocklicking.
function run_plotshocklicking_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotshocklicking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in lasso_classifier_pressmiss.
function lasso_classifier_pressmiss_Callback(hObject, eventdata, handles)
% hObject    handle to lasso_classifier_pressmiss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%lasso_classifier_pressmiss(handles.pathname,handles.filedate,handles.filename)
display('already calculated lasso');

% --- Executes on button press in run_plotlevertrials.
function run_plotlevertrials_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotlevertrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_plotlevertrials(handles.pathname,handles.filedate,handles.filename,handles.frametime);


% --- Executes on button press in run_plotshocktrials.
function run_plotshocktrials_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotshocktrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in run_plotPCA.
function run_plotPCA_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_plotPCA_v2(handles.pathname,handles.filedate,handles.filename,handles.frametime);


% --- Executes on button press in celltype_classifier_pressmiss.
function celltype_classifier_pressmiss_Callback(hObject, eventdata, handles)
% hObject    handle to celltype_classifier_pressmiss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%celltype_classifier_pressmiss(handles.pathname,handles.filedate,handles.filename,handles.frametime);
% run_plotLDA_5fold(handles.pathname,handles.filedate,handles.filename,handles.frametime);
celltype_classifier_pressmiss(handles.pathname,handles.filedate,handles.filename);

% --- Executes on button press in run_plotlevercells.
function run_plotlevercells_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotlevercells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
regression_all_neg(handles.pathname,handles.filedate,handles.filename,'lever');

% --- Executes on button press in run_plotrewardcells.
function run_plotrewardcells_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotrewardcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
regression_all_neg(handles.pathname,handles.filedate,handles.filename,'reward');

% --- Executes on button press in run_plotshockcells.
function run_plotshockcells_Callback(hObject, eventdata, handles)
% hObject    handle to run_plotshockcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
regression_all_neg(handles.pathname,handles.filedate,handles.filename,'shock');


% --- Executes on button press in run_regression_shockonlycells.
function run_regression_shockonlycells_Callback(hObject, eventdata, handles)
% hObject    handle to run_regression_shockonlycells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
regression_shockonlycells(handles.pathname,handles.filedate,handles.filename);

% --- Executes on button press in run_regression_barplot.
function run_regression_barplot_Callback(hObject, eventdata, handles)
% hObject    handle to run_regression_barplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
press_vs_miss_barplot(handles.pathname,handles.filedate,handles.filename);
