function varargout = koutu(varargin)
% KOUTU MATLAB code for koutu.fig
%      KOUTU, by itself, creates a new KOUTU or raises the existing
%      singleton*.
%
%      H = KOUTU returns the handle to a new KOUTU or the handle to
%      the existing singleton*.
%
%      KOUTU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KOUTU.M with the given input arguments.
%
%      KOUTU('Property','Value',...) creates a new KOUTU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before koutu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to koutu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help koutu

% Last Modified by GUIDE v2.5 12-Dec-2019 23:17:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @koutu_OpeningFcn, ...
                   'gui_OutputFcn',  @koutu_OutputFcn, ...
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


% --- Executes just before koutu is made visible.
function koutu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to koutu (see VARARGIN)

% Choose default command line output for koutu
handles.output = hObject;

global pos1 pos2 k m n_turn pic pic_size mask radius;
pos1=[];
pos2=[];
k=1000;
m=5;
n_turn=10;
pic = 1;

radius = 10;
axes(handles.axes1);
I = imread('.\data\1.jpg');
imshow(I);
pic_size = size(I);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes koutu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = koutu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% 画前景
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pos1 pic_size;
pos1 = [];
h1 = impoly(gca,'Closed',false);
foresub = getPosition(h1);
foresub = int32(foresub);
line(handles.axes1,foresub(:,1),foresub(:,2));
pos1 = sub2ind(pic_size,foresub(:,2),foresub(:,1));


%% 画背景
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pos2 pic_size;
pos2 = [];
h1 = impoly(gca,'Closed',false);
backsub = getPosition(h1);
backsub = int32(backsub);
line(handles.axes1,backsub(:,1),backsub(:,2),'color','red');
pos2 = sub2ind(pic_size,backsub(:,2),backsub(:,1));


%% 抠图
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pic pos1 pos2 k m n_turn mask;
if(pic == 1)
    I = imread('.\data\1.jpg');
elseif(pic == 2)
    I = imread('.\data\2.jpg');
else
    I = imread('.\data\3.jpg');
end
[L,N] = myslic(I, k, m, n_turn);
if(size(pos1,1)==0 || size(pos2,1)==0)
    return;
end
BW = lazysnapping(I,L,pos1,pos2);
mask = BW;
maskedImage = I;
maskedImage(repmat(~BW,[1 1 3])) = 0;
axes(handles.axes3);
imshow(maskedImage);

%% 选图片
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pic pic_size;
%cla(handles.axes1,'reset');
cla reset;
pic = get(handles.popupmenu1,'Value');
if(pic == 1)
    I = imread('.\data\1.jpg');
elseif(pic == 2)
    I = imread('.\data\2.jpg');
else
    I = imread('.\data\3.jpg');
end
pic_size = size(I);
axes(handles.axes1);
imshow(I)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% K个超像素
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k;
k = str2num(get(handles.edit1,'String'));
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% m值
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global m;
m = str2num(get(handles.edit2,'String'));
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% n_turn值
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global n_turn;
n_turn = str2num(get(handles.edit3,'String'));
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% 超像素过程显示 
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pic k m n_turn;
if(pic == 1)
    I = imread('.\data\1.jpg');
elseif(pic == 2)
    I = imread('.\data\2.jpg');
else
    I = imread('.\data\3.jpg');
end
%% 所用数据
[rows, cols, chan] = size(I);
img = rgb2lab(I);
S = sqrt(rows*cols/k);
S = ceil(S);
row_step = floor(rows/S);
col_step = floor(cols/S);
C = zeros(k,6);          % 1:3 mean Lab value; 4:5 x,y; 6 Num of Pixels
L = -ones(rows, cols);  
d = inf(rows, cols);  
% 初始化中心点
kk = 1;
for ii = 1:row_step
    for jj = 1:col_step
        rowc = round(S*(ii-0.5));
        colc = round(S*(jj-0.5));
        C(kk,1:3) = img(rowc,colc,:);
        C(kk,4) = rowc;
        C(kk,5) = colc;
        C(kk,6) = 0;
        kk = kk + 1;
    end
end
k = kk - 1;
% 开始迭代（一般不超过10次）
for n = 1:n_turn
    % assignment
    for kk = 1:k
        %搜索区域
        rmin = max(C(kk,4)-S, 1);   
        rmax = C(kk,4)+S; 
        if(rows-C(kk,4) < 2*S)
            rmax = rows;
        end
        cmin = max(C(kk,5)-S, 1);  
        cmax = C(kk,5)+S;
        if(cols-C(kk,5) < 2*S)
            cmax = cols;
        end
        
        for ii = rmin:rmax
            for jj =cmin:cmax
                dl = C(kk,1) - img(ii,jj,1);
                da = C(kk,2) - img(ii,jj,2);
                db = C(kk,3) - img(ii,jj,3);
                dx = C(kk,4) - ii;
                dy = C(kk,5) - jj;
                dc2 = dl^2 + da^2 + db^2;
                ds2 = dx^2 + dy^2;
                D = sqrt(dc2 + ds2 * m^2 / S^2);
                if(D < d(ii,jj))
                    d(ii,jj) = D;
                    L(ii,jj) = kk;
                end
            end
        end
    end
    %update
    C(:) = 0;
    for ii = 1:rows
         for jj = 1:cols
             C(L(ii,jj),1:5) = C(L(ii,jj),1:5) + [img(ii,jj,1) img(ii,jj,2) img(ii,jj,3) ii jj];
             C(L(ii,jj),6) = C(L(ii,jj),6) + 1;
         end
    end
    for kk = 1:k
        C(kk,1:5) = round(C(kk,1:5)/C(kk,6));
    end

    %显示中间过程图
    % 清除孤立像素
    for ii = 2:rows-1
        for jj = 2:cols-1
            this = L(ii,jj);
            same_num = (this==L(ii-1,jj-1)) + (this==L(ii-1,jj)) + (this==L(ii,jj-1)) + (this==L(ii+1 ,jj-1)) + (this==L(ii-1,jj+1)) + (this==L(ii,jj+1)) + (this==L(ii+1,jj)) + (this==L(ii+1,jj+1));
            if(same_num < 3.5)
                L(ii,jj) = L(ii-1,jj);
                C(this,6) = C(this,6) - 1;
            end
        end
    end
    % 重新整理L和C，得到N
    N = 0;
    to_zero = 0;%用来记kk前有几个超像素被清没了，L里标记减对应的数
    for kk = 1:k
      if(C(kk,6) ~= 0)
          N = N + 1;
      else
          to_zero = to_zero + 1;
      end
      C(kk,6) = to_zero;%功能变了，以前存包含个数，现在存之前空超像素个数
    end
    for ii = 1:rows
       for jj = 1:cols
           L(ii,jj) = L(ii,jj) - C(L(ii,jj),6);
       end
    end   
    BW = boundarymask(L);
    axes(handles.axes2);
    cla reset;
    imshow(imoverlay(I,BW,'cyan'),'InitialMagnification',67);
    hold on;
    scatter(C(:,5),C(:,4),'red','*');
    k = N;
end

%% 手动交互加前景
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask pic pic_size radius;
[col,row] = ginput;
[x,y] = meshgrid(1:pic_size(2),1:pic_size(1));
for i = 1:size(row,1)
    bw = sqrt((x-col(i,1)).^2 + (y-row(i,1)).^2) <= radius;
    mask = mask | bw;
end
if(pic == 1)
    I = imread('.\data\1.jpg');
elseif(pic == 2)
    I = imread('.\data\2.jpg');
else
    I = imread('.\data\3.jpg');
end
maskedImage = I;
maskedImage(repmat(~mask,[1 1 3])) = 0;
axes(handles.axes3);
imshow(maskedImage);


%% 手动交互减前景
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask pic pic_size radius;
[col,row] = ginput;
[x,y] = meshgrid(1:pic_size(2),1:pic_size(1));
for i = 1:size(row,1)
    bw = sqrt((x-col(i,1)).^2 + (y-row(i,1)).^2) <= radius;
    mask = mask & ~bw;
end
if(pic == 1)
    I = imread('.\data\1.jpg');
elseif(pic == 2)
    I = imread('.\data\2.jpg');
else
    I = imread('.\data\3.jpg');
end
maskedImage = I;
maskedImage(repmat(~mask,[1 1 3])) = 0;
axes(handles.axes3);
imshow(maskedImage);


%% 笔刷半径
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global radius;
radius = str2num(get(handles.edit4,'String'));
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
