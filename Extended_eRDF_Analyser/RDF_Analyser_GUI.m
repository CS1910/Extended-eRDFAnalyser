%% ------------------------------------------------------------------------
% ------------------------ INITIALISATION ---------------------------------
% -------------------------------------------------------------------------
function varargout = RDF_Analyser_GUI(varargin)
%
% RDF_ANALYSER_GUI MATLAB code for RDF_Analyser_GUI.fig
%
% RDF_ANALYSER_GUI is an interactive and integrated tool to extract reduced
% density functions (RDF) from electron diffraction patterns, and works
% for material compositions with upto 5 elements. For help on how to use
% the tool, please see the PDF User Manual. This program is free
% software, covered under the terms of GNU General Public License v3.
%
% -------------------------- Copyright (c) 2023 ---------------------------
% Christian Stenz
% 1st Institute of Physics (Ia)
% RWTH Aachen University
%
% Based on the source code from SHANMUGAM & BORISENKO (Copyright (c) 2017):
% Janaki Shanmugam & Konstantin B. Borisenko
% Electron Image Analysis Group, Department of Materials
% University of Oxford
% -------------------------------------------------------------------------
%
%      RDF_ANALYSER_GUI, by itself, creates a new RDF_ANALYSER_GUI or raises
%      the existing singleton*.
%
%      H = RDF_ANALYSER_GUI returns the handle to a new RDF_ANALYSER_GUI or the
%      handle to the existing singleton*.
%
%      RDF_ANALYSER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RDF_ANALYSER_GUI.M with the given input arguments.
%
%      RDF_ANALYSER_GUI('Property','Value',...) creates a new RDF_ANALYSER_GUI or
%      raises the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RDF_Analyser_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RDF_Analyser_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RDF_Analyser_GUI

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RDF_Analyser_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @RDF_Analyser_GUI_OutputFcn, ...
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

% -------------------- SET (DEFAULT) PARAMETERS ---------------------------
function RDF_Analyser_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

set(0,'defaulttextinterpreter','latex')

% This function has no output args, see OutputFcn.
% varargin   command line arguments to RDF_Analyser_GUI (see VARARGIN)
dbstop if error
% Choose default command line output for RDF_Analyser_GUI
handles.output = hObject;

% --- Default values to be used
set(handles.Tab1,'Value',1); % depressed Tab1
% in case edit box callback functions are not executed

assignin('base', 'SelectedFitCurve',1);
handles.Firsttime = 1;  % sets to 0 when "RDF Plot" Tab is opened for the first time and function is plotted
handles.mx = 1;
% default checkboxes
handles.load_calib = 1;
handles.overwrite  = 0;
handles.save       = 0;

handles.ds = 0.00224;
set(handles.edit_ds, 'String', num2str(handles.ds))

handles.elem1 = 33; %Ge
handles.elem2 = 1; %none
handles.elem3 = 1; %none
handles.elem4 = 1; %none
handles.elem5 = 1; %none
handles.e1 = 1;
handles.e2 = 0;
handles.e3 = 0;
handles.e4 = 0;
handles.e5 = 0;
handles.q_fix = 0;
handles.dq = 0.1;
handles.N = 10;
handles.dN = 1;
handles.edit_dN = 1;
handles.edit_damping = 0.5;
handles.paramK = load('Kirkland_2010.txt','-ascii');
handles.param_val = 2;
handles.rmax = 10;
handles.Perc = 30;

set(handles.Tab2,'Value',1); % depressed Tab2
% radio button groups - default selection
set(handles.uibuttongroup_DP,'SelectedObject',handles.Amorphous); % centre optimisation routine

% display default values used for fitting
set(handles.text_q_fit, 'String', handles.q_fix);
set(handles.text_N, 'String', handles.N);
set(handles.text_damping, 'String', handles.edit_damping);
set(handles.popup_param, 'Value', 2); %Kirkland
set(handles.Element1, 'Value', 34); %Ge
set(handles.Element2, 'Value', 2); %none
set(handles.Element3, 'Value', 2); %none
set(handles.Element4, 'Value', 2); %none
set(handles.Element5, 'Value', 2); %none

% Update handles structure
guidata(hObject, handles);
% -------------------------------------------------------------------------
function varargout = RDF_Analyser_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
varargout{1} = handles.output;
% -------------------------------------------------------------------------
function About_ClickedCallback(hObject, eventdata, handles)
(msgbox...
    ({'eRDF Analyser is distributed in the hope that it will be useful,'...
    'but without any warranty. This program is free software and'...
    'you are welcome to redistribute it under certain conditions.'...
    'See the GNU General Public License for more details:'...
    '(http://www.gnu.org/licenses/).'...
    ''...
    'Copyright (c) 2016'...
    'J Shanmugam, KB Borisenko'},'About eRDF Analyser'));
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% --------------------------- TAB CONTROL ---------------------------------
% -------------------------------------------------------------------------
function Tab1_Callback(hObject, eventdata, handles)
set(handles.Panel1,'Visible','On');
set(handles.Panel2,'Visible','Off');
set(handles.Tab2,'Value',0);
function Tab2_Callback(hObject, eventdata, handles)
set(handles.Panel1,'Visible','Off');
set(handles.Panel2,'Visible','On');
if isfield(handles, 'fileNumber')
    foldername = handles.datfname;
    CurveString = char(foldername{:});
    set(handles.popSelectCurve, 'String', CurveString);
    set(handles.popSelectCurve, 'Value', 1);
    Button_Iq_Callback(handles.Button_Iq, eventdata, handles);
end
set(handles.Tab1,'Value',0);



%% ------------------------------------------------------------------------
% ---------------------------- CALIBRATION --------------------------------
% -------------------------------------------------------------------------
function text_result_eds_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '');
function text_result_ds_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '');
function checkbox_load_calibration_Callback(hObject, eventdata, handles)
handles.load_calib = get(hObject,'Value');
guidata(hObject,handles)

% -------------------- CALIBRATION VIA 2D SAED PATTERN --------------------
function Button_Calibration_Data_Callback(hObject, eventdata, handles)
% Choose data from input file (diffraction pattern image as text file)
[fname,pname] = uigetfile({'*.tif','TIF';'*.*','All files (*.*)';},...
    'Choose input data file','Test_data\GST225',...
    'MultiSelect', 'off');

if iscell(fname) ~= 1
    if fname == 0
        % User clicked the Cancel button.
        return;
    end
    fname = char2cell(fname);
end

% mx needed for the equidistant Colors in the plots.
if size(fname,2) == 1
    mx0 = 1;
else
    mx0 = size(fname,2) - 1;
end

fileNumber = size(fname,2);
for fNo = 1:fileNumber

    fname0 = [pname,fname{fNo}];
    parts = strsplit(fname0, '\');
    FileTag = [parts{end-2}, '\..\', parts{end}(1:end-4)];
    addpath(pname);
    rehash toolboxcache;

    % check if already calculated
    Res_ds_Name = sprintf('%s%s_Result_ds.txt', pname, fname{fNo}(1:end-4));
    if isfile(Res_ds_Name) && handles.load_calib

        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.DataLines = [2, Inf];
        opts.Delimiter = "\t";
        opts.VariableNames = ["dsA1px", "uncertaintyA1px"];
        opts.VariableTypes = ["double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        results = readtable(Res_ds_Name, opts);
        results = table2array(results);
        clear opts

        ds  = results(1,1);
        eds_sys = results(1,2);

        guidata(hObject,handles)
        handles.ds = ds;
        set(handles.edit_ds, 'String', num2str(ds))
        set(handles.text_result_ds, 'String', ['(' num2str(ds)])
        set(handles.text_result_eds, 'String', [num2str(eds_sys) ')'])
        set(handles.text_calibration_file, 'String', FileTag)
        guidata(hObject,handles)
        return
    else

        [pathstr,name,ext] = fileparts(fname0); %#ok<ASGLU>
        guidata(hObject,handles)

        % Load 8bit RGB image (series 1)-> '1'
        dp = imread([pname,fname{fNo}], 1);
        % RBG layer is again divided into 3 layers for R/G/B, respectively.
        % In black/white image only one layer needs to be loaded to obtain all
        % information, since all layers are identical.
        dp = dp(:,:,1);  % use layer 1 ("R")
        dp = double(dp);

        % Load 16bit greyscale image (series 2)-> '2'
        dpraw16b = imread([pname,fname{fNo}], 2);
        dpraw16b = double(dpraw16b);

        % If no 16bit layer is available, use the 8bit layer instead:
%         dpraw16b = dp;

        nx = size(dp,1);  % rows
        ny = size(dp,2);  % coloumns

        % For better visualisation:
        % Apply median filter to remove single pixel noise
        dp   = medfilt2(dp);
        dp16 = medfilt2(dpraw16b);
        ig   = log(abs(dp)+1);
        ig16 = log(abs(dpraw16b)+1);
        plim = mat2gray(ig);

        % Show SAED Pattern (8bit layer)
        DPfig  = figure('Name',['Diffraction Pattern ', fname{fNo}],...
            'NumberTitle','off');
        imshow(imadjust(plim));

        % ---------------------------------------------------------------------
        %% Masking parts of the image
        % ---------------------------------------------------------------------

        % Create mask
        TotalMask = zeros(size(dp));

        %% ------------ Automated masking ---------------
%         % Mask scale bar
%         scaleMask = zeros(size(dp));
%         scaleMask(end-80:end,1:370) = ones(81,370);  % lower left corner
%         scaleMask(end-80:end,end-275:end) = ones(81,276);  % lower right corner
%         scaleMask = logical(scaleMask);
%         TotalMask = TotalMask|scaleMask(size(dp));

        
        % Mask beam stop (automatically)
        beamstopMask = importdata('MaskShape.txt');  % import beam stop mask shape
        DPMaskSectionX = (nx/2-150:nx/2+150);  % nx rows
        DPMaskSectionY = (ny/2-170:ny/2+170);  % ny coloumns  -> central part around intensity maximum

        im = imbinarize(dp(DPMaskSectionX,DPMaskSectionY),15);  % create true b/w image of beamstop in image in selected area (15 = threshold)
        mask = beamstopMask(DPMaskSectionX,DPMaskSectionY);  % crop same area from the mask shape matrix

        % compute cross-correlation of both matrices and find maximum
        crr = xcorr2(-1.*im+1,mask);
        [~,indexm] = max(crr(:));

        % convert index of maximum to x- & y-indices and compute the actual
        % shift of the mask in order to mask the beam stop
        crrxdim = size(crr,2);
        crrydim = size(crr,1);
        yshift = ceil(indexm./crrydim);
        yshift = yshift-size(im,2);
        xshift = mod(indexm,crrydim);
        if xshift == 0
            xshift = crrydim;
        end
        xshift = xshift-size(im,1);

        % shift the mask and correct the overlapping edges
        beamstopMask = logical(circshift(beamstopMask,[xshift yshift]));  % (row,coloumn)
        if yshift > 0
            beamstopMask(:,1:yshift) = 0;
        elseif yshift < 0
            beamstopMask(972+xshift:1030+xshift,end+yshift:end) = 1;
        end

        % Burn mask into image by setting it to NaN (0) wherever the mask is true.
        TotalMask = TotalMask|beamstopMask;
        dpraw16b(TotalMask) = NaN;
        dp16(TotalMask)     = NaN;
        ig16(TotalMask)     = 0;
        dp(TotalMask)       = NaN;
        ig(TotalMask)       = 0;

        % Display the masked image (16bit layer)
        plim16 = mat2gray(ig16);
        figure(DPfig);
        imshow(imadjust(plim16));


        %% ------------ Add handdrawn masks ---------------


%         % RECTANGLE
%         % Alert user to mask beam stop with RECTANGLES
%         uiwait(msgbox({'Click and drag to draw a rectangle to mask beamstop';
%             'or any desired region. If not, double-click anywhere in figure.'},...
%             'Masking'));
%         FreeMask = 1; %--- while loop (to create additional Masks)
%         while FreeMask == 1
% 
%             % Create rectangular ROI
%             hFH = drawrectangle();
%             binaryMask_rect = hFH.createMask();
% 
%             TotalMask = TotalMask|binaryMask_rect;
% 
%             % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%             dpraw16b(TotalMask) = NaN;
%             dp16(TotalMask)     = NaN;
%             ig16(TotalMask)     = 0;
%             dp(TotalMask)       = NaN;
%             ig(TotalMask)       = 0;
% 
%             % Display the masked image (16bit layer)
%             plim16 = mat2gray(ig16);
%             figure(DPfig);
%             imshow(imadjust(plim16));
% 
%             % save shape of mask
%             filename0 = [pname, 'Mask.txt'];
%             fid=fopen(filename0,'wt');
%             for ii = 1:size(TotalMask,1)
%                 fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%                 fprintf(fid,'\n');
%             end
% 
%             % ask user for additional masks
%             AddMask = questdlg(['Do you want to mask additional regions ' ...
%                 'in rectangular shape?'],...
%                 'Additional Masking',...
%                 'No','Yes',...
%                 'No');
%             switch AddMask
%                 case 'No'
%                     FreeMask = 0;
%                     continue %--- exit while loop
%                 case 'Yes'
%                     FreeMask = 1;
%             end
%         end
% 
% 
% 
%         % ELLIPSE
%         % Alert user to mask beam stop with ELLIPSE
%         uiwait(msgbox({'Click and drag to draw an ellipse to mask beamstop';
%             'or any desired region. If not, double-click anywhere in figure.'},...
%             'Masking'));
%         FreeMask = 1; %--- while loop (to create additional Masks)
%         while FreeMask == 1
% 
%             % Create rectangular ROI
%             hFH = drawellipse();
%             binaryMask_elli = hFH.createMask();
% 
%             TotalMask = TotalMask|binaryMask_elli;
% 
%             % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%             dpraw16b(TotalMask) = NaN;
%             dp16(TotalMask)     = NaN;
%             ig16(TotalMask)     = 0;
%             dp(TotalMask)       = NaN;
%             ig(TotalMask)       = 0;
% 
%             % Display the masked image (16bit layer)
%             plim16 = mat2gray(ig16);
%             figure(DPfig);
%             imshow(imadjust(plim16));
% 
%             % save shape of mask
%             filename0 = [pname, 'Mask.txt'];
%             fid=fopen(filename0,'wt');
%             for ii = 1:size(TotalMask,1)
%                 fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%                 fprintf(fid,'\n');
%             end
% 
%             % ask user for additional masks
%             AddMask = questdlg(['Do you want to mask additional regions ' ...
%                 'in ellipsoid shape?'],...
%                 'Additional Masking',...
%                 'No','Yes',...
%                 'No');
%             switch AddMask
%                 case 'No'
%                     FreeMask = 0;
%                     continue %--- exit while loop
%                 case 'Yes'
%                     FreeMask = 1;
%             end
%         end
% 
% 
% 
%         % POLYGON
%         % Alert user to mask beam stop with POLYGONS
%         uiwait(msgbox({'Click and draw a polygon to mask beamstop';
%             'or any desired region. If not, double-click anywhere in figure.'},...
%             'Masking'));
%         FreeMask = 1; %--- while loop (to create additional Masks)
%         while FreeMask == 1
% 
%             % Create rectangular ROI
%             hFH = drawpolygon();
%             binaryMask_poly = hFH.createMask();
% 
%             TotalMask = TotalMask|binaryMask_poly;
% 
%             % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%             dpraw16b(TotalMask) = NaN;
%             dp16(TotalMask)     = NaN;
%             ig16(TotalMask)     = 0;
%             dp(TotalMask)       = NaN;
%             ig(TotalMask)       = 0;
% 
%             % Display the masked image (16bit layer)
%             plim16 = mat2gray(ig16);
%             figure(DPfig);
%             imshow(imadjust(plim16));
% 
%             % save shape of mask
%             filename0 = [pname, 'Mask.txt'];
%             fid=fopen(filename0,'wt');
%             for ii = 1:size(TotalMask,1)
%                 fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%                 fprintf(fid,'\n');
%             end
% 
%             % ask user for additional masks
%             AddMask = questdlg(['Do you want to mask additional regions ' ...
%                 'in polygon shape?'],...
%                 'Additional Masking',...
%                 'No','Yes',...
%                 'No');
%             switch AddMask
%                 case 'No'
%                     FreeMask = 0;
%                     continue %--- exit while loop
%                 case 'Yes'
%                     FreeMask = 1;
%             end
%         end




        % -------------------------------------------------------
        %% Find centre
        % -------------------------------------------------------
        
        % Parameters for automatic circe identification
        binary_threshold = 2.0;  % try 0.1 to 2.0 or more
        min_circle_dia   = 300;  % minimum circle diameter in pixels that can be identified as central circular object
        max_circle_dia   = 2000;  % maximum circle diameter in pixels that can be identified as central circular object
        OptGridSize      = 6;  % 13x13 pixel grid that is scanned for optimal center position around the central circular object
        
        OPT = 1; %--- while optimisation loop
        while OPT == 1
            %% Automated center finder

            % find first estimate of center
            im = dp16;
            Img_max		= max(max(im));
            Img_min		= min(min(im));
            Img_Delta	= Img_max - Img_min;
            Img_Avg_pre	= (im - Img_min)./Img_Delta;

            % Smooth image
            Img_Avg_pre2 = wiener2(Img_Avg_pre, [2 2]);
            Img_Avg      = imgaussfilt(Img_Avg_pre2 ,2);

            % Binarize image
            Img_level	= graythresh(Img_Avg);
            Img_bin	    = imbinarize(Img_Avg, Img_level*binary_threshold);
            imshow(Img_bin);


            % Possibly working to fill the diffraction rings with the missing
            % beam stop part to make an identification of the interrupted ring
            % pattern as a circle in the next step more probable. To try
            % uncomment the following lines.
            %             for ir = 1:ny/2
            %                 i = ny-ir;
            %                 if sum(Img_bin(:,i)) > 4 && sum(Img_bin(:,i-1)) > 5
            %                     j=1;
            %                     while Img_bin(j,i) == 0 && j < nx
            %                         j=j+1;
            %                     end
            %                     while Img_bin(j,i) == 1 && j < nx
            %                         j=j+1;
            %                     end
            %                     countZeros = 0;
            %                     jz = j;
            %                     while Img_bin(jz,i) == 0 && jz < nx
            %                         countZeros = countZeros + 1;
            %                         jz = jz + 1;
            %                     end
            %                     if countZeros >= 22 && countZeros <= 75
            %                         Img_bin(j:j+countZeros,i-400:i) = 1;
            %                         break;
            %                     end
            %                 end
            %             end
            %
            %             fillCircleMask = beamstopMask;
            %             fillCircleMask(:,i:end) = logical(0);
            %             Img_bin = Img_bin | fillCircleMask;
            %             imshow(Img_bin);


            % Find circles / convex hull
            CH_objects = bwconvhull(Img_bin,'objects');
            Img_bin = Img_bin | CH_objects;
            imshow(Img_bin);

            stats = regionprops('table', Img_bin,...
                'Centroid', 'Eccentricity', 'EquivDiameter');
            hBright = viscircles(stats.Centroid(:,:),...
                stats.EquivDiameter(:,:)/2,...
                'Color','b');

            % Filter found circles to diameters between 300px & 2000px
            stats = stats(stats.EquivDiameter > min_circle_dia ...
                        & stats.EquivDiameter < max_circle_dia, :);
            delete(hBright);
            hBright = viscircles(stats.Centroid(:,:),...
                stats.EquivDiameter(:,:)/2, 'Color','r');

            xc = stats.Centroid(1,1);
            yc = stats.Centroid(1,2);
            diameter = stats.EquivDiameter(1,1);

            radius = 0.5*diameter;
            xMin = xc - radius;
            yMin = yc - radius;

            % create ellipse
            hEllipse = imellipse(gca,[xMin, yMin, diameter, diameter]);
            hEllipse.setFixedAspectRatioMode('true');
            pos = hEllipse.getPosition;
            delete(hEllipse)

            clear stats xc yc diameter radius xMin yMin


            %% User input of center
% 
%             % Plot DP (16 bit)
%             figure(DPfig);
%             clear gcf
%             displaydp = dpraw16b;
%             displaydp(TotalMask)= min(min(displaydp))*0.8;
%             displaydp = mat2gray(displaydp);
%             imshow(imadjust(displaydp));
%             colormap(gca,colorcube);
%             
%             % create ellipse
%             hEllipse = imellipse(gca,[nx/2, ny/2, nx/3, ny/3]);
%             hEllipse.setFixedAspectRatioMode( 'true' );
% 
%             % Alert user to adjust ellipse position and double click when finished
%             h = helpdlg({'Move and resize marker to fit one of the inner contours',...
%                 '(in black / blue / green). Try to be as accurate as possible.','',...
%                 'Double-click inside ellipse once finished.'},...
%                 'Adjust and double-click marker');
% 
%             % wait for double-click -> get position
%             wait(hEllipse);
%             if ishandle(DPfig) == 0 %user closes DPfig, terminate and return to main GUI
%                 return
%             end
%             pos = hEllipse.getPosition;
% 
%             % close helpdlg if user doesn't
%             try close (h)
%             catch
%             end
            %% -----------------------------------------------------------

            
            % Plot the center and ellipse
            figure(DPfig);
            hold on;
            rad = 0.5*pos(3);
            dm = pos(3);
            xcentre = pos(1)+rad;
            ycentre = pos(2)+rad;
            Circle_User = rectangle( ...
                'Position',[xcentre-rad, ycentre-rad, dm, dm], ...
                'Curvature', [1,1], 'EdgeColor', 'white', 'LineWidth', 2);
            Centre_User = plot(xcentre, ycentre, 'r+', 'LineWidth', 1,...
                'MarkerSize', 20, 'HandleVisibility', 'off');
            delete(hEllipse);

            % case 'Polycrystalline'
            % Plot DP (16 bit)
            figure(DPfig);
            clear gcf
            displaydp = dpraw16b;
            displaydp(TotalMask)= min(min(displaydp))*0.8;
            displaydp = mat2gray(displaydp);
            imshow(imadjust(displaydp));
            colormap(gca,colorcube);

            prompt = {'Enter number of projections:',...
                'Enter distance of contour from edge (in pixels)',...
                'Enter size of grid scan (in pixels)'};
            % default values
            def = {'100','200','15'};
            opt_input = inputdlg(prompt,'Optimisation parameters',1,def);
            if isempty(opt_input)
                % User clicked cancel
                OPT = 0; % exit optimisation routine
                figure(DPfig);
                delete(Circle_User);
                delete(Centre_User);
                continue
            end
            % input values into variables
            nnp      = str2double(opt_input{1});
            dedge    = str2double(opt_input{2});
            maxshift = str2double(opt_input{3});
            figure(DPfig);
            delete(Circle_User);
            % ----------------------------------------------------
            % Optimize the circle position globally with upscaled precision
            % ----------------------------------------------------
            scale = 10;

            % Redefine initial outline with more points between selected angles
            amin=-85*pi/180;
            amax=85*pi/180;
            alpha=amin:(amax-amin)/(nnp-1):amax;
            outline=zeros(nnp,2);
            outline(:,1)=rad*cos(alpha)+xcentre;
            outline(:,2)=rad*sin(alpha)+ycentre;

            % Coordinates of the line scans
            sl=zeros(nnp,1);
            cc=zeros(nnp,1);

            for nn=1:nnp
                vv=outline(nn,1)-xcentre;
                sl(nn)=(outline(nn,2)-ycentre)/vv;
                cc(nn)=scale*(outline(nn,2)-sl(nn)*outline(nn,1));
            end

            % The largest full circle in original scaling
            % Distances from centre close the dedge distance from the edge

            edist=zeros(4,1);
            edist(1)=abs(nx-xcentre)-dedge;
            edist(2)=abs(ny-ycentre)-dedge;
            edist(3)=xcentre-dedge;
            edist(4)=ycentre-dedge;
            p2=min(edist);
            Circle_p2 = rectangle('Position',[xcentre-p2, ycentre-p2, 2*p2, 2*p2], ...
                'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
            rad_opt = rad; % for plotting of Circle_User2

            xo1=zeros(nnp,1);
            xo2=zeros(nnp,1);
            yo1=zeros(nnp,1);
            yo2=zeros(nnp,1);
            xi1=zeros(nnp,1);
            xi2=zeros(nnp,1);
            yi1=zeros(nnp,1);
            yi2=zeros(nnp,1);

            % Sections in upscaled version
            for nn=1:nnp
                % Inner circle
                aqi=sl(nn)*sl(nn)+1;
                bqi=2*(sl(nn)*cc(nn)-sl(nn)*ycentre*scale-xcentre*scale);
                cqi=xcentre*xcentre*scale*scale+cc(nn)*cc(nn)-2*cc(nn)*ycentre*scale+...
                    ycentre*ycentre*scale*scale-rad*rad*scale*scale;
                ddi=bqi*bqi-4*aqi*cqi;
                xi1(nn)=round((-bqi+sqrt(ddi))/(2*aqi));
                xi2(nn)=round((-bqi-sqrt(ddi))/(2*aqi));
                yi1(nn)=round(sl(nn)*xi1(nn)+cc(nn));
                yi2(nn)=round(sl(nn)*xi2(nn)+cc(nn));

                % Outer circle
                aqo=sl(nn)*sl(nn)+1;
                bqo=2*(sl(nn)*cc(nn)-sl(nn)*ycentre*scale-xcentre*scale);
                cqo=xcentre*xcentre*scale*scale+cc(nn)*cc(nn)-2*cc(nn)*ycentre*scale+...
                    ycentre*ycentre*scale*scale-p2*p2*scale*scale;
                ddo=bqo*bqo-4*aqo*cqo;
                xo1(nn)=round((-bqo+sqrt(ddo))/(2*aqo));
                xo2(nn)=round((-bqo-sqrt(ddo))/(2*aqo));
                yo1(nn)=round(sl(nn)*xo1(nn)+cc(nn));
                yo2(nn)=round(sl(nn)*xo2(nn)+cc(nn));
            end
            % -------------------------------------------------------
            % Optimise centre using centrosymmetric line profiles
            % -------------------------------------------------------
            % Upscaled calculations
            nsample=1000*scale;

            lineprofile1=zeros(nsample,nnp);
            lineprofile2=zeros(nsample,nnp);

            if(maxshift > dedge)
                maxshift=dedge-1;
            end
            ssum=zeros(2*maxshift+1,2*maxshift+1);

            % initial positions
            xi1_o=xi1;
            xo1_o=xo1;
            xi2_o=xi2;
            xo2_o=xo2;
            cc_o=cc;

            % Initialise waitbar
            ProgBar = waitbar(0,'Please wait...','Name','Optimising centre',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(ProgBar,'canceling',0)
            full = length(-maxshift:maxshift);
            count = 0;
            br = 0; % flag for canceling

            for ii=-maxshift:maxshift
                count = count+1;
                for jj=-maxshift:maxshift
                    xi1=xi1_o+jj*scale;
                    xo1=xo1_o+jj*scale;
                    xi2=xi2_o+jj*scale;
                    xo2=xo2_o+jj*scale;
                    cc=cc_o+ii*scale-jj*sl*scale;

                    for nn=1:nnp
                        if (xi1(nn) < xo1(nn))
                            xx=xi1(nn):(xo1(nn)-xi1(nn))/(nsample-1):xo1(nn);
                        else
                            xx=xo1(nn):(xi1(nn)-xo1(nn))/(nsample-1):xi1(nn);
                        end
                        yy=sl(nn)*xx+cc(nn);
                        % Downscale and get values
                        xxd=round(xx/scale);
                        yyd=round(yy/scale);

                        indx=sub2ind(size(dp),yyd,xxd);
                        lineprofile1(:,nn)=dp(indx);

                        if (xi2(nn) < xo2(nn))
                            xx=xi2(nn):(xo2(nn)-xi2(nn))/(nsample-1):xo2(nn);
                        else
                            xx=xo2(nn):(xi2(nn)-xo2(nn))/(nsample-1):xi2(nn);
                        end
                        yy=sl(nn)*xx+cc(nn);
                        % Downscale and get values
                        xxd=round(xx/scale);
                        yyd=round(yy/scale);

                        indx=sub2ind(size(dp),yyd,xxd);
                        sindx=max(size(indx));
                        indxt=zeros(1,nsample);
                        indxt(1:1:sindx)=indx(sindx:-1:1);
                        lineprofile2(:,nn)=dp(indxt);

                        diff=lineprofile1(:,nn)-lineprofile2(:,nn);
                        if(sum(isnan(diff)) == 0)
                            ssum(ii+maxshift+1,jj+maxshift+1)=ssum(ii+maxshift+1,jj+maxshift+1)+sum(diff.*diff);
                        end
                    end
                end
                % Check for Cancel button press
                if getappdata(ProgBar,'canceling')
                    br = 1;
                    break;
                end
                % Update waitbar
                waitbar(count/full,ProgBar,'Please wait...');
            end
            delete(ProgBar);
            if br == 1
                figure(DPfig);
                delete(Centre_User)
                delete(Circle_p2)
                OPT = 1;
                continue
            end

            % centre shift
            [~,optxys]=min(ssum(:));    % [optval,optxys]
            [optxs,optys]=ind2sub(size(ssum),optxys);
            optxshift=optys-maxshift-1;
            optyshift=optxs-maxshift-1;

            xcentre_opt=xcentre+optxshift;
            ycentre_opt=ycentre+optyshift;

            % -------------------------------------------------------
            % Plot the sum of optimised line profiles to test the fit
            % -------------------------------------------------------
            % Optimised line profile positions
            xi1=xi1_o+optxshift*scale;
            xo1=xo1_o+optxshift*scale;
            xi2=xi2_o+optxshift*scale;
            xo2=xo2_o+optxshift*scale;
            cc=cc_o+optyshift*scale-optxshift*sl*scale;

            sum_lineprofile1=zeros(nsample,1);
            sum_lineprofile2=zeros(nsample,1);

            for nn=1:nnp
                if (xi1(nn) < xo1(nn))
                    xx=xi1(nn):(xo1(nn)-xi1(nn))/(nsample-1):xo1(nn);
                else
                    xx=xo1(nn):(xi1(nn)-xo1(nn))/(nsample-1):xi1(nn);
                end
                yy=sl(nn)*xx+cc(nn);
                % Downscale and get values
                xxd=round(xx/scale);
                yyd=round(yy/scale);

                indx=sub2ind(size(dp),yyd,xxd);
                lineprofile1(:,nn)=dp(indx);

                if (xi2(nn) < xo2(nn))
                    xx=xi2(nn):(xo2(nn)-xi2(nn))/(nsample-1):xo2(nn);
                else
                    xx=xo2(nn):(xi2(nn)-xo2(nn))/(nsample-1):xi2(nn);
                end
                yy=sl(nn)*xx+cc(nn);
                % Downscale and get values
                xxd=round(xx/scale);
                yyd=round(yy/scale);

                indx=sub2ind(size(dp),yyd,xxd);
                sindx=max(size(indx));
                indxt=zeros(1,nsample);
                indxt(1:1:sindx)=indx(sindx:-1:1);
                lineprofile2(:,nn)=dp(indxt);

                if(sum(isnan(lineprofile1(:,nn))+isnan(lineprofile2(:,nn))) == 0)
                    sum_lineprofile1=sum_lineprofile1+lineprofile1(:,nn);
                    sum_lineprofile2=sum_lineprofile2+lineprofile2(:,nn);
                end
            end

            figure('Name','Centrefinder fit','NumberTitle','off');
            plot(sum_lineprofile1);
            hold all;
            plot(sum_lineprofile2);
            plot(sum_lineprofile1-sum_lineprofile2);
            legend('User line profile','Optimised line profile','Difference');

            if (abs(xcentre_opt-xcentre) == maxshift || abs(ycentre_opt-ycentre) == maxshift)
                %                 uiwait(msgbox({'Optimisation not successful!';
                %                     'Increase maxshift or provide better initial guess.'},...
                %                     'Error','error'));
                OPT = 1;

                % Show image (16bit)
                figure(DPfig);
                delete(Centre_User);
                delete(Circle_p2);
                OptGridSize = OptGridSize+2;
                try
                    delete(Circle_p1);
                catch
                end
                continue
            end

            % Plot DP (16 bit)
            figure(DPfig);
            clear gcf
            displaydp = dpraw16b;
            displaydp(TotalMask)= min(min(displaydp))*0.8;
            displaydp = mat2gray(displaydp);
            imshow(imadjust(displaydp));

            % Plot the initial circle
            Circle_User2 = rectangle(...
                'Position',[xcentre-rad, ycentre-rad, 2*rad, 2*rad], ...
                'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
            % Plot the optimised circle and centre
            Centre_Opt = plot(xcentre_opt, ycentre_opt, 'g+', 'LineWidth', 1,...
                'MarkerSize', 20, 'HandleVisibility', 'off');
            Circle_Opt = rectangle(...
                'Position',[xcentre_opt-rad_opt, ycentre_opt-rad_opt, 2*rad_opt, 2*rad_opt], ...
                'Curvature', [1,1], 'EdgeColor', 'g', 'LineWidth', 2);
            try
                Circle_p1_Opt = rectangle(...
                    'Position',[xcentre_opt-p1, ycentre_opt-p1, 2*p1, 2*p1], ...
                    'Curvature', [1,1], 'EdgeColor', 'm', 'LineWidth', 2);
            catch
            end
            Circle_p2_Opt = rectangle(...
                'Position',[xcentre_opt-p2, ycentre_opt-p2, 2*p2, 2*p2], ...
                'Curvature', [1,1], 'EdgeColor', 'm', 'LineWidth', 2);

            guidata(hObject,handles)

            xcentre = xcentre_opt;
            ycentre = ycentre_opt;
            OPT = 0; %--- end while (optimisation) loop

            % delete rings
            figure(DPfig);
            delete(Centre_User);
            delete(Centre_Opt);
            delete(Circle_User2);
            delete(Circle_Opt);
            delete(Circle_p2);
            delete(Circle_p2_Opt);
        end
        % -------------------------------------------------------
        %% Calculate azimuthal average and variance with centre (xcentre, ycentre)
        % -------------------------------------------------------

        % Distances from corners to the centre of the diffraction pattern
        
        % Original scale
        cdist=zeros(4,1);
        cdist(1)=sqrt((nx-xcentre)^2+(ny-ycentre)^2);
        cdist(2)=sqrt((xcentre)^2+(ny-ycentre)^2);
        cdist(3)=sqrt((nx-xcentre)^2+(ycentre)^2);
        cdist(4)=sqrt((xcentre)^2+(ycentre)^2);

        azsize=ceil(max(cdist));
        azav=zeros(azsize,1);
        nazav=zeros(azsize,1);
        m2=zeros(azsize,1);
        delta=zeros(azsize,1);

        % Mean & variance
        for xx=1:nx
            for yy=1:ny
                if (~isnan(dpraw16b(yy,xx)))
                    kk=ceil(sqrt((xx-xcentre)^2+(yy-ycentre)^2));
                    nazav(kk)=nazav(kk)+1;
                    delta(kk)=dpraw16b(yy,xx)-azav(kk);
                    azav(kk)=azav(kk)+delta(kk)./nazav(kk);
                    m2(kk)=m2(kk)+delta(kk)*(dpraw16b(yy,xx)-azav(kk));
                end
            end
        end
        azvar=m2./nazav;
        % Normalized variance -------
        % nazvar=azvar./(azav.*azav);
        guidata(hObject,handles)
        % -------------------------------------------------------
        % Plot and save average and variance data
        % -------------------------------------------------------
        % x-axis for plots (pixel index)
        pix_xax = linspace(1,azsize,azsize);
        q_xax = pix_xax';

        % Plot DP (16 bit)
        figure(DPfig);
        clear gcf
        displaydp = dpraw16b;
        displaydp(TotalMask)= min(min(displaydp))*0.8;
        displaydp = mat2gray(displaydp);
        imshow(imadjust(displaydp));
        colormap(gca,colorcube);

        rectangle('Position',[xcentre 1 nx-xcentre ycentre-1],...
            'FaceColor',[1 1 1 0.9]);

        plot(pix_xax+xcentre,ycentre-1-azav*(0.9*ycentre)/max(azav),...
            'Color', Color(3), 'DisplayName','azimuthal average');
        plot(pix_xax+xcentre,ycentre-1-azvar*(0.9*ycentre)/max(azvar),...
            'Color', Color(11), 'DisplayName','azimuthal variance',...
            'LineStyle', '--');
        legend('show');

        % Azimuthal Average
        if fNo == 1
            AzimAv = figure('Name','Azimuthal Average',...
                'NumberTitle','off');
        else
            figure(AzimAv);
            hold on
        end
        plot(q_xax,azav, 'DisplayName', fname{fNo}, 'Color',...
            Color(3+7.5/mx0*(fNo-1)));
        xlabel('Pixel');
        ylabel('Intensity');
        legend('show', 'Interpreter', 'none');

        % Azimuthal Variance
        if fNo == 1
            AzimVar = figure('Name','Azimuthal Variance',...
                'NumberTitle','off');
        else
            figure(AzimVar);
            hold on
        end
        plot(q_xax,azvar, 'DisplayName', fname{fNo}, 'Color',...
            Color(3+7.5/mx0*(fNo-1)));
        xlabel('Pixel');
        ylabel('Variance');
        legend('show', 'Interpreter', 'none');

        folder = pname;
        addpath(folder);
        guidata(hObject,handles)

        azav_name = sprintf('%s/%s_azav.txt',folder,name);
        save (azav_name,'azav','-ASCII');
        azvar_name = sprintf('%s/%s_azvar.txt',folder,name);
        save (azvar_name,'azvar','-ASCII');

        % -------------------------------------------------------
        % compute average of all input images
        if fNo == 1
            sum_azav_name = ['_', fname{fNo}(1:4)];
            sum_azav = azav;
        else
            sum_azav_name = [sum_azav_name, '_', fname{fNo}(1:4)];
            lenDif = size(azav,1) - size(sum_azav,1);
            if lenDif > 0
                sum_azav = [sum_azav; zeros(lenDif,1)] + azav;
            elseif lenDif < 0
                sum_azav = [azav; zeros(abs(lenDif),1)] + sum_azav;
            else
                sum_azav = sum_azav + azav;
            end
        end

        delete(Circle_User)
        delete(Circle_Opt)
        delete(Circle_p2_Opt)
        saveas(DPfig, [pname , '\', fname{fNo}(1:end-4),'_DPCenter.fig'])
    end
end
sum_azav = sum_azav./fileNumber;
sum_azav_Name = sprintf(...
    ['%s_sum_azav', sum_azav_name,'.txt'],folder);
save(sum_azav_Name,'sum_azav','-ASCII');

AzimAv_Name = sprintf(...
    ['%s_AzimAv', sum_azav_name,'.fig'],folder);
AzimVar_Name = sprintf(...
    ['%s_AzimVar', sum_azav_name,'.fig'],folder);


%% Calibrate peaks to d-spacings

sum_azav_ct = sum_azav(65:end);  % avoid peak from intensity maximum at center

% find peaks
[pks,locs,w,prom] = findpeaks(sum_azav_ct);
[pks,locs,w,prom] = findpeaks(sum_azav_ct, 'MinPeakProminence', 0.005*max(prom));

locs = locs + 65;

% sort by peak heigth
[pks_sort, indx] = sort(pks, 'descend');
locs_sort = locs(indx);
w_sort    = w(indx);
prom_sort = prom(indx);

% only use 5 highest peaks
pks_sort  = pks_sort(1:5);
locs_sort = locs_sort(1:5);
w_sort    = w_sort(1:5);
prom_sort = prom_sort(1:5);

% sort by pixel value
[px_locs, indx] = sort(locs_sort);
pks_sort   = pks_sort(indx);
w_sort     = w_sort(indx);
prom_sort  = prom_sort(indx);

figure(AzimAv);
hold on
plot(px_locs,pks_sort,'x', 'Color', 'r', 'DisplayName', 'Fitted Peaks')
saveas(AzimAv, AzimAv_Name)
saveas(AzimVar, AzimVar_Name)


% read in d-spacing file;  values in Angstrom
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["dhkl", "uncertainty", "hkl"];
opts.VariableTypes = ["double", "double", "string"];
opts.ConsecutiveDelimitersRule = "join";
opts = setvaropts(opts, "hkl", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "hkl", "EmptyFieldRule", "auto");
Au_T = readtable("Au_mean_std.txt", opts);
clear opts

dhkl_fit = Au_T.dhkl([1,2,3,4,7]);
dhkl_err = Au_T.uncertainty([1,2,3,4,7]);
hkl      = Au_T.hkl([1,2,3,4,7]);

% reciprocal lattice parameter
Khkl     = 1./dhkl_fit;
Khkl_err = dhkl_err./(dhkl_fit).^2;

% fitting with errors --> systematic error eds_sys
CalPlot = figure('Name', 'Calibration Fit', 'NumberTitle', 'off');
[P,SP] = linfitxy(px_locs, Khkl, 1., Khkl_err, 'Verbosity', 0);
grid on
xlabel('Pixel');
ylabel('$K_{hkl}$ (\AA$$^{-1}$$)', 'Interpreter', 'Latex');

% save fit
CalPlot_Name = sprintf('%s%s_Calibration_Fit.fig', folder, name);
saveas(CalPlot, CalPlot_Name);

ds      = P(1);
eds_sys = SP(1);

guidata(hObject,handles)
handles.ds = ds;
set(handles.edit_ds, 'String', num2str(ds))
set(handles.text_result_ds, 'String', ['(' num2str(ds)])
set(handles.text_result_eds, 'String', [num2str(eds_sys) ')'])
set(handles.text_calibration_file, 'String', FileTag)
guidata(hObject,handles)

%% Write file with results for loading a calibration factor
Res_ds = strings(2,2);
Res_ds(1,1) = 'ds (A^(-1)/px)';
Res_ds(1,2) = 'uncertainty (A^(-1)/px)';
Res_ds(2,1) = num2str(ds);
Res_ds(2,2) = num2str(eds_sys);

Res_ds_Name = sprintf('%s%s_Result_ds.txt', folder, name);
fid = fopen(Res_ds_Name,'wt');
for row = 1:size(Res_ds,1)
    fprintf(fid, repmat('%s\t',1,size(Res_ds,2)-1), Res_ds{row,1:end-1});
    fprintf(fid, '%s \n', Res_ds{row,end});
end
fclose(fid);

% ------------------------ CALIBRATION FACTOR -----------------------------
function edit_ds_Callback(hObject, eventdata, handles)
ds = str2double(get(hObject,'String'));
if isnan(ds)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
handles.ds = ds;
guidata(hObject,handles)
function edit_ds_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', num2str(0.00224))
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% ------------------------------------------------------------------------
% ------------------------ OPEN 2D SAED TIF -------------------------------
% -------------------------------------------------------------------------
function button_OpenDP_Callback(hObject, eventdata, handles)

% Choose data from input file (diffraction pattern image as text file)
[fname,pname] = uigetfile({'*.tif','TIF';'*.*','All files (*.*)';},...
    'Choose input data file','Test_data\GST225',...
    'MultiSelect', 'on');

if iscell(fname) ~= 1
    if fname == 0
        % User clicked the Cancel button.
        return;
    end
    fname = char2cell(fname);
end

% mx needed for the equidistant Colors in the plots.
if size(fname,2) == 1
    mx0 = 1;
else
    mx0 = size(fname,2) - 1;
end

fileNumber = size(fname,2);
for fNo = 1:fileNumber

    fname0 = [pname,fname{fNo}];
    addpath(pname);
    rehash toolboxcache;
    handles.pname = pname;

    [pathstr,name,ext] = fileparts(fname0); %#ok<ASGLU>
    handles.DPfname = name;
    guidata(hObject,handles)


    % TIF contains 2 layers (in our case): 8bit RGB and 16bit greyscale.
    % 8bit image is better for contrast/brightness manipulation and
    % identifiction of the beam blank.
    % In contrast to 16bit image, the 8bit image also includes the scale bar.
    % Identification of the center, azimuthal averaging and calculation of
    % the intensity profile is done using the 16bit greyscale, as it
    % contains much more information.
    
    % Load 8bit RGB image (series 1)-> '1'
    dp = imread([pname,fname{fNo}], 1);
    % RBG layer is again divided into 3 layers for R/G/B, respectively.
    % In black/white image only one layer needs to be loaded to obtain all 
    % information, since all layers are identical.
    dp = dp(:,:,1);  % use layer 1 ("R")
    dp = double(dp);
    

    % Load 16bit greyscale image (series 2)-> '2'
    dpraw16b = imread([pname,fname{fNo}], 2);
    dpraw16b = double(dpraw16b);

    % If no 16bit layer is available, use the 8bit layer instead:
%     dpraw16b = dp;

    nx = size(dp,1);  % rows
    ny = size(dp,2);  % coloumns

    % For better visualisation:
    % Apply median filter to remove single pixel noise
    dp   = medfilt2(dp);
    dp16 = medfilt2(dpraw16b);
    ig   = log(abs(dp)+1);
    ig16 = log(abs(dpraw16b)+1);
    plim = mat2gray(ig);


    % Show SAED Pattern (8bit layer)
    DPfig  = figure('Name',['Diffraction Pattern ', fname{fNo}],...
        'NumberTitle','off');
    imshow(imadjust(plim));

    % ---------------------------------------------------------------------
    %% Masking parts of the image
    % ---------------------------------------------------------------------

    % Create mask
    TotalMask = zeros(size(dp));

    %% ------------ Automated masking ---------------
%     % Mask scale bar
%     scaleMask = zeros(size(dp));
%     scaleMask(end-80:end,1:370) = ones(81,370);  % lower left corner
%     scaleMask(end-80:end,end-275:end) = ones(81,276);  % lower right corner
%     scaleMask = logical(scaleMask);
%     TotalMask = TotalMask|scaleMask(size(dp));

    % Mask beam stop (automtically)

    % -----------------------------------------------------------------
    beamstopMask = importdata('MaskShape.txt');  % import beam stop mask shape
    DPMaskSectionX = (nx/2-150:nx/2+150);  % nx rows
    DPMaskSectionY = (ny/2-170:ny/2+170);  % ny coloumns  -> central part around intensity maximum

    im = imbinarize(dp(DPMaskSectionX,DPMaskSectionY),80);  % create true b/w image of beamstop in image in selected area
    mask = beamstopMask(DPMaskSectionX,DPMaskSectionY);  % crop same area from the mask shape matrix

    % compute cross-correlation of both matrices and find maximum
    crr = xcorr2(-1.*im+1,mask);
    [~,indexm] = max(crr(:));

    % convert index of maximum to x- & y-indices and compute the actual
    % shift of the mask in order to mask the beam stop
    crrxdim = size(crr,2);
    crrydim = size(crr,1);
    yshift = ceil(indexm./crrydim);
    yshift = yshift-size(im,2);
    xshift = mod(indexm,crrydim);
    if xshift == 0
        xshift = crrydim;
    end
    xshift = xshift-size(im,1);

    % shift the mask and correct the overlapping edges
    beamstopMask = logical(circshift(beamstopMask,[xshift yshift]));  % (row,coloumn)
    if yshift > 0
        beamstopMask(:,1:yshift) = 0;
    elseif yshift < 0
        beamstopMask(972+xshift:1030+xshift,end+yshift:end) = 1;
    end

    % Burn mask into image by setting it to NaN (0) wherever the mask is true.
    TotalMask = TotalMask|beamstopMask;
    dpraw16b(TotalMask) = NaN;
    dp16(TotalMask)     = NaN;
    ig16(TotalMask)     = 0;
    dp(TotalMask)       = NaN;
    ig(TotalMask)       = 0;

    % Display the masked image (16bit layer)
    plim16 = mat2gray(ig16);
    figure(DPfig);
    imshow(imadjust(plim16));


    %% ------------ Add handdrawn masks ---------------

 
%     % RECTANGLE
%     % Alert user to mask beam stop with RECTANGLES
%     uiwait(msgbox({'Click and drag to draw a rectangle to mask beamstop';
%         'or any desired region. If not, double-click anywhere in figure.'},...
%         'Masking'));
%     FreeMask = 1; %--- while loop (to create additional Masks)
%     while FreeMask == 1
% 
%         % Create rectangular ROI
%         hFH = drawrectangle();
%         binaryMask_rect = hFH.createMask();
% 
%         TotalMask = TotalMask|binaryMask_rect;
% 
%         % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%         dpraw16b(TotalMask) = NaN;
%         dp16(TotalMask)     = NaN;
%         ig16(TotalMask)     = 0;
%         dp(TotalMask)       = NaN;
%         ig(TotalMask)       = 0;
% 
%         % Display the masked image (16bit layer)
%         plim16 = mat2gray(ig16);
%         figure(DPfig);
%         imshow(imadjust(plim16));
% 
%         % save shape of mask
%         filename0 = [pname, 'Mask.txt'];
%         fid=fopen(filename0,'wt');
%         for ii = 1:size(TotalMask,1)
%             fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%             fprintf(fid,'\n');
%         end
% 
%         % ask user for additional masks
%         AddMask = questdlg(['Do you want to mask additional regions ' ...
%             'in rectangular shape?'],...
%             'Additional Masking',...
%             'No','Yes',...
%             'No');
%         switch AddMask
%             case 'No'
%                 FreeMask = 0;
%                 continue %--- exit while loop
%             case 'Yes'
%                 FreeMask = 1;
%         end
%     end
% 
% 
% 
%     % ELLIPSE
%     % Alert user to mask beam stop with ELLIPSE
%     uiwait(msgbox({'Click and drag to draw an ellipse to mask beamstop';
%         'or any desired region. If not, double-click anywhere in figure.'},...
%         'Masking'));
%     FreeMask = 1; %--- while loop (to create additional Masks)
%     while FreeMask == 1
% 
%         % Create rectangular ROI
%         hFH = drawellipse();
%         binaryMask_elli = hFH.createMask();
% 
%         TotalMask = TotalMask|binaryMask_elli;
% 
%         % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%         dpraw16b(TotalMask) = NaN;
%         dp16(TotalMask)     = NaN;
%         ig16(TotalMask)     = 0;
%         dp(TotalMask)       = NaN;
%         ig(TotalMask)       = 0;
% 
%         % Display the masked image (16bit layer)
%         plim16 = mat2gray(ig16);
%         figure(DPfig);
%         imshow(imadjust(plim16));
% 
%         % save shape of mask
%         filename0 = [pname, 'Mask.txt'];
%         fid=fopen(filename0,'wt');
%         for ii = 1:size(TotalMask,1)
%             fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%             fprintf(fid,'\n');
%         end
% 
%         % ask user for additional masks
%         AddMask = questdlg(['Do you want to mask additional regions ' ...
%             'in ellipsoid shape?'],...
%             'Additional Masking',...
%             'No','Yes',...
%             'No');
%         switch AddMask
%             case 'No'
%                 FreeMask = 0;
%                 continue %--- exit while loop
%             case 'Yes'
%                 FreeMask = 1;
%         end
%     end
% 
% 
% 
%     % POLYGON
%     % Alert user to mask beam stop with POLYGONS
%     uiwait(msgbox({'Click and draw a polygon to mask beamstop';
%         'or any desired region. If not, double-click anywhere in figure.'},...
%         'Masking'));
%     FreeMask = 1; %--- while loop (to create additional Masks)
%     while FreeMask == 1
% 
%         % Create rectangular ROI
%         hFH = drawpolygon();
%         binaryMask_poly = hFH.createMask();
% 
%         TotalMask = TotalMask|binaryMask_poly;
% 
%         % Burn mask into image by setting it to NaN (0) wherever the mask is true.
%         dpraw16b(TotalMask) = NaN;
%         dp16(TotalMask)     = NaN;
%         ig16(TotalMask)     = 0;
%         dp(TotalMask)       = NaN;
%         ig(TotalMask)       = 0;
% 
%         % Display the masked image (16bit layer)
%         plim16 = mat2gray(ig16);
%         figure(DPfig);
%         imshow(imadjust(plim16));
% 
%         % save shape of mask
%         filename0 = [pname, 'Mask.txt'];
%         fid=fopen(filename0,'wt');
%         for ii = 1:size(TotalMask,1)
%             fprintf(fid,'%6.1f \t',TotalMask(ii,:));
%             fprintf(fid,'\n');
%         end
% 
%         % ask user for additional masks
%         AddMask = questdlg(['Do you want to mask additional regions ' ...
%             'in polygon shape?'],...
%             'Additional Masking',...
%             'No','Yes',...
%             'No');
%         switch AddMask
%             case 'No'
%                 FreeMask = 0;
%                 continue %--- exit while loop
%             case 'Yes'
%                 FreeMask = 1;
%         end
%     end
%     


    % -------------------------------------------------------
    %% Find centre
    % -------------------------------------------------------

    % Parameters for circe identification
    binary_threshold = 0.5;  % try 0.1 to 2.0 or more
    min_circle_dia = 300;  % minimum circle diameter in pixels that can be identified as central circular object
    max_circle_dia = 2000;  % maximum circle diameter in pixels that can be identified as central circular object
    OptGridSize = 6;  % 13x13 pixel grid that is scanned for optimal center position around the central circular object

    OPT = 1; %--- while optimisation loop
    while OPT == 1
        %% Automated center finder

        % find first estimate of center
        im = dp16;
        Img_max		= max(max(im));
        Img_min		= min(min(im));
        Img_Delta	= Img_max - Img_min;
        Img_Avg_pre	= (im - Img_min)./Img_Delta;

        % Smooth image
        Img_Avg_pre2  = wiener2(Img_Avg_pre, [2 2]);
        Img_Avg       = imgaussfilt(Img_Avg_pre2 ,2);

        % Binarize image
        Img_level	= graythresh(Img_Avg);
        Img_bin		= imbinarize(Img_Avg, Img_level*binary_threshold);
        imshow(Img_bin);


        % Possibly working to fill the diffraction rings with the missing
        % beam stop part to make an identification of the interrupted ring
        % pattern as a circle in the next step more probable. To try
        % uncomment the following lines.
        % for ir = 1:ny/2
        %     i = ny-ir;
        %     if sum(Img_bin(:,i)) > 4 && sum(Img_bin(:,i-1)) > 5
        %         j=1;
        %         while Img_bin(j,i) == 0 && j < nx
        %             j=j+1;
        %         end
        %         while Img_bin(j,i) == 1 && j < nx
        %             j=j+1;
        %        end
        %         countZeros = 0;
        %         jz = j;
        %         while Img_bin(jz,i) == 0 && jz < nx
        %             countZeros = countZeros + 1;
        %             jz = jz + 1;
        %         end
        %         if countZeros >= 22 && countZeros <= 75
        %             Img_bin(j:j+countZeros,i-400:i) = 1;
        %             break;
        %         end
        %     end
        % end
        %
        % fillCircleMask = beamstopMask;
        % fillCircleMask(:,i:end) = logical(0);
        % Img_bin = Img_bin | fillCircleMask;
        % imshow(Img_bin);

        % Find circles / convex hull
        CH_objects = bwconvhull(Img_bin,'objects');
        Img_bin = Img_bin | CH_objects;
        imshow(Img_bin);

        % Find circles
        stats = regionprops('table', Img_bin, ...
            'Centroid', 'Eccentricity', 'EquivDiameter');
        hBright = viscircles(stats.Centroid(:,:), ...
            stats.EquivDiameter(:,:)/2,'Color','b');

        % Filter found circles to diameters between 300px & 2000px
        stats = stats(stats.EquivDiameter > min_circle_dia ...
                    & stats.EquivDiameter < max_circle_dia, :);
        delete(hBright);
        hBright = viscircles(stats.Centroid(:,:), ...
            stats.EquivDiameter(:,:)/2,'Color','r');

        xc = stats.Centroid(1,1);
        yc = stats.Centroid(1,2);
        diameter = stats.EquivDiameter(1,1);

        radius = 0.5*diameter;
        xMin = xc - radius;
        yMin = yc - radius;
        
        % create ellipse
        hEllipse = imellipse(gca,[xMin, yMin, diameter, diameter]);
        hEllipse.setFixedAspectRatioMode( 'true' );
        pos = hEllipse.getPosition;
        delete(hEllipse)

        clear stats xc yc diameter radius xMin yMin

        
        %% User input of center
% 
%         % Plot DP (16 bit)
%         figure(DPfig);
%         clear gcf
%         displaydp = dpraw16b;
%         displaydp(TotalMask)= min(min(displaydp))*0.8;
%         displaydp = mat2gray(displaydp);
%         imshow(imadjust(displaydp));
%         colormap(gca,colorcube);
% 
%         % create ellipse
%         hEllipse = imellipse(gca,[nx/2, ny/2, nx/3, ny/3]);
%         hEllipse.setFixedAspectRatioMode( 'true' );
% 
%         % Alert user to adjust ellipse position and double click when finished
%         h = helpdlg({'Move and resize marker to fit one of the inner contours',...
%             '(in black / blue / green). Try to be as accurate as possible.','',...
%             'Double-click inside ellipse once finished.'},...
%             'Adjust and double-click marker');
% 
%         % wait for double-click -> get position
%         wait(hEllipse);
%         if ishandle(DPfig) == 0 %user closes DPfig, terminate and return to main GUI
%             return
%         end
%         pos = hEllipse.getPosition;
% 
%         % close helpdlg if user doesn't
%         try close (h)
%         catch
%         end
        %% ---------------------------------------------------------------

        % Plot the center and ellipse
        figure(DPfig);
        hold on;
        rad = 0.5*pos(3);
        dm = pos(3);
        xcentre = pos(1)+rad;
        ycentre = pos(2)+rad;
        Circle_User = rectangle( ...
            'Position',[xcentre-rad, ycentre-rad, dm, dm], ...
            'Curvature', [1,1], 'EdgeColor', 'white', 'LineWidth', 2);
        Centre_User = plot(xcentre, ycentre, 'r+', 'LineWidth', 1,...
            'MarkerSize', 20, 'HandleVisibility', 'off');
        delete(hEllipse);

        DPtype = get(get(handles.uibuttongroup_DP,'SelectedObject'),'Tag');
        switch DPtype
            case 'Polycrystalline' %----------------------------------
                % CentreOpt routine #1: radial lines
                % ----------------------------------------------------
                % Input dialog for optimisation parameters
                prompt = {'Enter number of projections:',...
                    'Enter distance of contour from edge (in pixels)',...
                    'Enter size of grid scan (in pixels)'};
                % default values
                def = {'100','200','15'};
                opt_input = inputdlg(prompt,'Optimisation parameters',1,def);
                if isempty(opt_input)
                    % User clicked cancel
                    OPT = 0; % exit optimisation routine
                    figure(DPfig);
                    delete(Circle_User);
                    delete(Centre_User);
                    continue
                end
                % input values into variables
                nnp      = str2double(opt_input{1});
                dedge    = str2double(opt_input{2});
                maxshift = str2double(opt_input{3});
                figure(DPfig);
                delete(Circle_User);
                % ----------------------------------------------------
                % Optimize the circle position globally with upscaled precision
                % ----------------------------------------------------
                scale = 10;

                % Redefine initial outline with more points between selected angles
                amin=-85*pi/180;
                amax=85*pi/180;
                alpha=amin:(amax-amin)/(nnp-1):amax;
                outline=zeros(nnp,2);
                outline(:,1)=rad*cos(alpha)+xcentre;
                outline(:,2)=rad*sin(alpha)+ycentre;

                % Coordinates of the line scans
                sl=zeros(nnp,1);
                cc=zeros(nnp,1);

                for nn=1:nnp
                    vv=outline(nn,1)-xcentre;
                    sl(nn)=(outline(nn,2)-ycentre)/vv;
                    cc(nn)=scale*(outline(nn,2)-sl(nn)*outline(nn,1));
                end

                % The largest full circle in original scaling
                % Distances from centre close the dedge distance from the edge

                edist=zeros(4,1);
                edist(1)=abs(nx-xcentre)-dedge;
                edist(2)=abs(ny-ycentre)-dedge;
                edist(3)=xcentre-dedge;
                edist(4)=ycentre-dedge;
                p2=min(edist);
                Circle_p2 = rectangle('Position',[xcentre-p2, ycentre-p2, 2*p2, 2*p2], ...
                    'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
                rad_opt = rad; % for plotting of Circle_User2

                xo1=zeros(nnp,1);
                xo2=zeros(nnp,1);
                yo1=zeros(nnp,1);
                yo2=zeros(nnp,1);
                xi1=zeros(nnp,1);
                xi2=zeros(nnp,1);
                yi1=zeros(nnp,1);
                yi2=zeros(nnp,1);

                % Sections in upscaled version
                for nn=1:nnp
                    % Inner circle
                    aqi=sl(nn)*sl(nn)+1;
                    bqi=2*(sl(nn)*cc(nn)-sl(nn)*ycentre*scale-xcentre*scale);
                    cqi=xcentre*xcentre*scale*scale+cc(nn)*cc(nn)-2*cc(nn)*ycentre*scale+...
                        ycentre*ycentre*scale*scale-rad*rad*scale*scale;
                    ddi=bqi*bqi-4*aqi*cqi;
                    xi1(nn)=round((-bqi+sqrt(ddi))/(2*aqi));
                    xi2(nn)=round((-bqi-sqrt(ddi))/(2*aqi));
                    yi1(nn)=round(sl(nn)*xi1(nn)+cc(nn));
                    yi2(nn)=round(sl(nn)*xi2(nn)+cc(nn));

                    % Outer circle
                    aqo=sl(nn)*sl(nn)+1;
                    bqo=2*(sl(nn)*cc(nn)-sl(nn)*ycentre*scale-xcentre*scale);
                    cqo=xcentre*xcentre*scale*scale+cc(nn)*cc(nn)-2*cc(nn)*ycentre*scale+...
                        ycentre*ycentre*scale*scale-p2*p2*scale*scale;
                    ddo=bqo*bqo-4*aqo*cqo;
                    xo1(nn)=round((-bqo+sqrt(ddo))/(2*aqo));
                    xo2(nn)=round((-bqo-sqrt(ddo))/(2*aqo));
                    yo1(nn)=round(sl(nn)*xo1(nn)+cc(nn));
                    yo2(nn)=round(sl(nn)*xo2(nn)+cc(nn));
                end
                % -------------------------------------------------------
                % Optimise centre using centrosymmetric line profiles
                % -------------------------------------------------------
                % Upscaled calculations
                nsample=1000*scale;

                lineprofile1=zeros(nsample,nnp);
                lineprofile2=zeros(nsample,nnp);

                if(maxshift > dedge)
                    maxshift=dedge-1;
                end
                ssum=zeros(2*maxshift+1,2*maxshift+1);

                % initial positions
                xi1_o=xi1;
                xo1_o=xo1;
                xi2_o=xi2;
                xo2_o=xo2;
                cc_o=cc;

                % Initialise waitbar
                ProgBar = waitbar(0,'Please wait...','Name','Optimising centre',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(ProgBar,'canceling',0)
                full = length(-maxshift:maxshift);
                count = 0;
                br = 0; % flag for canceling

                for ii=-maxshift:maxshift
                    count = count+1;
                    for jj=-maxshift:maxshift
                        xi1=xi1_o+jj*scale;
                        xo1=xo1_o+jj*scale;
                        xi2=xi2_o+jj*scale;
                        xo2=xo2_o+jj*scale;
                        cc=cc_o+ii*scale-jj*sl*scale;

                        for nn=1:nnp
                            if (xi1(nn) < xo1(nn))
                                xx=xi1(nn):(xo1(nn)-xi1(nn))/(nsample-1):xo1(nn);
                            else
                                xx=xo1(nn):(xi1(nn)-xo1(nn))/(nsample-1):xi1(nn);
                            end
                            yy=sl(nn)*xx+cc(nn);
                            % Downscale and get values
                            xxd=round(xx/scale);
                            yyd=round(yy/scale);

                            indx=sub2ind(size(dp),yyd,xxd);
                            lineprofile1(:,nn)=dp(indx);

                            if (xi2(nn) < xo2(nn))
                                xx=xi2(nn):(xo2(nn)-xi2(nn))/(nsample-1):xo2(nn);
                            else
                                xx=xo2(nn):(xi2(nn)-xo2(nn))/(nsample-1):xi2(nn);
                            end
                            yy=sl(nn)*xx+cc(nn);
                            % Downscale and get values
                            xxd=round(xx/scale);
                            yyd=round(yy/scale);

                            indx=sub2ind(size(dp),yyd,xxd);
                            sindx=max(size(indx));
                            indxt=zeros(1,nsample);
                            indxt(1:1:sindx)=indx(sindx:-1:1);
                            lineprofile2(:,nn)=dp(indxt);

                            diff=lineprofile1(:,nn)-lineprofile2(:,nn);
                            if(sum(isnan(diff)) == 0)
                                ssum(ii+maxshift+1,jj+maxshift+1)=ssum(ii+maxshift+1,jj+maxshift+1)+sum(diff.*diff);
                            end
                        end
                    end
                    % Check for Cancel button press
                    if getappdata(ProgBar,'canceling')
                        br = 1;
                        break;
                    end
                    % Update waitbar
                    waitbar(count/full,ProgBar,'Please wait...');
                end
                delete(ProgBar);
                if br == 1
                    figure(DPfig);
                    delete(Centre_User)
                    delete(Circle_p2)
                    OPT = 1;
                    continue
                end

                % centre shift
                [~,optxys]=min(ssum(:));    % [optval,optxys]
                [optxs,optys]=ind2sub(size(ssum),optxys);
                optxshift=optys-maxshift-1;
                optyshift=optxs-maxshift-1;

                xcentre_opt=xcentre+optxshift;
                ycentre_opt=ycentre+optyshift;

                % -------------------------------------------------------
                % Plot the sum of optimised line profiles to test the fit
                % -------------------------------------------------------
                % Optimised line profile positions
                xi1=xi1_o+optxshift*scale;
                xo1=xo1_o+optxshift*scale;
                xi2=xi2_o+optxshift*scale;
                xo2=xo2_o+optxshift*scale;
                cc=cc_o+optyshift*scale-optxshift*sl*scale;

                sum_lineprofile1=zeros(nsample,1);
                sum_lineprofile2=zeros(nsample,1);

                for nn=1:nnp
                    if (xi1(nn) < xo1(nn))
                        xx=xi1(nn):(xo1(nn)-xi1(nn))/(nsample-1):xo1(nn);
                    else
                        xx=xo1(nn):(xi1(nn)-xo1(nn))/(nsample-1):xi1(nn);
                    end
                    yy=sl(nn)*xx+cc(nn);
                    % Downscale and get values
                    xxd=round(xx/scale);
                    yyd=round(yy/scale);

                    indx=sub2ind(size(dp),yyd,xxd);
                    lineprofile1(:,nn)=dp(indx);

                    if (xi2(nn) < xo2(nn))
                        xx=xi2(nn):(xo2(nn)-xi2(nn))/(nsample-1):xo2(nn);
                    else
                        xx=xo2(nn):(xi2(nn)-xo2(nn))/(nsample-1):xi2(nn);
                    end
                    yy=sl(nn)*xx+cc(nn);
                    % Downscale and get values
                    xxd=round(xx/scale);
                    yyd=round(yy/scale);

                    indx=sub2ind(size(dp),yyd,xxd);
                    sindx=max(size(indx));
                    indxt=zeros(1,nsample);
                    indxt(1:1:sindx)=indx(sindx:-1:1);
                    lineprofile2(:,nn)=dp(indxt);

                    if(sum(isnan(lineprofile1(:,nn))+isnan(lineprofile2(:,nn))) == 0)
                        sum_lineprofile1=sum_lineprofile1+lineprofile1(:,nn);
                        sum_lineprofile2=sum_lineprofile2+lineprofile2(:,nn);
                    end
                end

                figure('Name','Centrefinder fit','NumberTitle','off');
                plot(sum_lineprofile1);
                hold all;
                plot(sum_lineprofile2);
                plot(sum_lineprofile1-sum_lineprofile2);
                legend('User line profile','Optimised line profile','Difference');


            case 'Amorphous'  %---------------------------------------
                % CentreOpt routine #2: by minimising sum of azimuthal variance
                % ----------------------------------------------------
                % Input dialog for optimisation parameters
                %                     prompt = {'Enter size of grid scan (in pixels)'};
                %                     def = {'10'}; % default value
                %                     opt_input = inputdlg(prompt,'Optimisation parameters',1,def);
                %                     if isempty(opt_input)
                %                         % User clicked cancel
                %                         OPT = 0; % exit optimisation routine
                %                         figure(DPfig);
                %                         delete(Circle_User);
                %                         delete(Centre_User);
                %                         continue
                %                     end

                opt_input = {num2str(OptGridSize)};

                % input values into variables
                maxshift = str2double(opt_input{1});
                figure(DPfig);
                delete(Circle_User);
                % ----------------------------------------------------
                % Define limits for centre optimisation
                p1 = ceil(rad*0.75);
                p2 = ceil(rad+p1);
                Circle_p1 = rectangle('Position',[xcentre-p1, ycentre-p1, 2*p1, 2*p1], ...
                    'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
                Circle_p2 = rectangle('Position',[xcentre-p2, ycentre-p2, 2*p2, 2*p2], ...
                    'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
                % ----------------------------------------------------
                % Initialise waitbar
                ProgBar = waitbar(0,'Please wait...','Name','Optimising centre',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(ProgBar,'canceling',0);
                full = length(-maxshift:maxshift);
                count = 0;
                br = 0; % flag for canceling

                ssum=zeros(2*maxshift+1,2*maxshift+1);
                azavsize=p2-p1+1;

                for ii=-maxshift:maxshift
                    count = count + 1;
                    for jj=-maxshift:maxshift
                        azav=zeros(azavsize,1);
                        nazav=zeros(azavsize,1);
                        m2=zeros(azavsize,1);
                        delta=zeros(azavsize,1);
                        % Mean & variance
                        for xx=1:nx
                            for yy=1:ny
                                if (~isnan(dpraw16b(xx,yy)))
                                    kk=ceil(sqrt((xx-ycentre+ii)^2+(yy-xcentre+jj)^2));
                                    if (kk >= p1 && kk <= p2)
                                        pp=kk-p1+1;
                                        nazav(pp)=nazav(pp)+1;
                                        delta(pp)=dpraw16b(xx,yy)-azav(pp);
                                        azav(pp)=azav(pp)+delta(pp)./nazav(pp);
                                        m2(pp)=m2(pp)+delta(pp)*(dpraw16b(xx,yy)-azav(pp));
                                    end
                                end
                            end
                        end
                        % Check for Cancel button press in waitbar
                        if getappdata(ProgBar,'canceling')
                            br = 1;
                            break;
                        end
                        azvar=m2./nazav;
                        icount=-maxshift:maxshift==ii;
                        jcount=-maxshift:maxshift==jj;
                        ssum(icount,jcount)=sum(azvar.^2);
                    end
                    if br == 1
                        break;
                    end
                    % Update waitbar
                    waitbar(count/full,ProgBar,'Please wait...');
                end
                delete(ProgBar);
                if br == 1
                    figure(DPfig);
                    delete(Centre_User)
                    delete(Circle_p1)
                    delete(Circle_p2)
                    OPT = 1;
                    continue
                end

                [optval,optxys]=min(ssum(:)); %#ok<ASGLU>

                [optxs,optys]=ind2sub(size(ssum),optxys);
                optxshift=optys-maxshift-1;
                optyshift=optxs-maxshift-1;

                xcentre_opt=xcentre-optxshift;
                ycentre_opt=ycentre-optyshift;
                rad_opt = rad;

        end

        if (abs(xcentre_opt-xcentre) == maxshift || abs(ycentre_opt-ycentre) == maxshift)
            %                 uiwait(msgbox({'Optimisation not successful!';
            %                     'Increase maxshift or provide better initial guess.'},...
            %                     'Error','error'));
            OPT = 1;

            % Show image (16bit)
            figure(DPfig);
            delete(Centre_User);
            delete(Circle_p2);
            OptGridSize = OptGridSize+2;
            try
                delete(Circle_p1);
            catch
            end
            continue
        end

        % Plot DP (16 bit)
        figure(DPfig);
        clear gcf
        displaydp = dpraw16b;
        displaydp(TotalMask)= min(min(displaydp))*0.8;
        displaydp = mat2gray(displaydp);
        imshow(imadjust(displaydp));

        % Plot the initial circle
        Circle_User2 = rectangle(...
            'Position',[xcentre-rad, ycentre-rad, 2*rad, 2*rad], ...
            'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 2);
        % Plot the optimised circle and centre
        Centre_Opt = plot(xcentre_opt, ycentre_opt, 'g+', 'LineWidth', 1,...
            'MarkerSize', 20, 'HandleVisibility', 'off');
        Circle_Opt = rectangle(...
            'Position',[xcentre_opt-rad_opt, ycentre_opt-rad_opt, 2*rad_opt, 2*rad_opt], ...
            'Curvature', [1,1], 'EdgeColor', 'g', 'LineWidth', 2);
        try
            Circle_p1_Opt = rectangle(...
                'Position',[xcentre_opt-p1, ycentre_opt-p1, 2*p1, 2*p1], ...
                'Curvature', [1,1], 'EdgeColor', 'm', 'LineWidth', 2);
        catch
        end
        Circle_p2_Opt = rectangle(...
            'Position',[xcentre_opt-p2, ycentre_opt-p2, 2*p2, 2*p2], ...
            'Curvature', [1,1], 'EdgeColor', 'm', 'LineWidth', 2);

        guidata(hObject,handles)

        xcentre = xcentre_opt;
        ycentre = ycentre_opt;
        OPT = 0; %--- end while (optimisation) loop
        
        % delete rings
        figure(DPfig);
        delete(Centre_User);
        delete(Centre_Opt);
        delete(Circle_User2);
        delete(Circle_Opt);
        delete(Circle_p2);
        delete(Circle_p2_Opt);
    end
    % -------------------------------------------------------
    %% Calculate azimuthal average and variance with centre (xcentre, ycentre)
    % -------------------------------------------------------
    % Distances from corners to the centre of the diffraction pattern
    % Original scale

    cdist=zeros(4,1);
    cdist(1)=sqrt((nx-xcentre)^2+(ny-ycentre)^2);
    cdist(2)=sqrt((xcentre)^2+(ny-ycentre)^2);
    cdist(3)=sqrt((nx-xcentre)^2+(ycentre)^2);
    cdist(4)=sqrt((xcentre)^2+(ycentre)^2);

    azsize=ceil(max(cdist));
    azav=zeros(azsize,1);
    nazav=zeros(azsize,1);
    m2=zeros(azsize,1);
    delta=zeros(azsize,1);

    % Mean & variance
    for xx=1:nx
        for yy=1:ny
            if (~isnan(dpraw16b(yy,xx)))
                kk=ceil(sqrt((xx-xcentre)^2+(yy-ycentre)^2));
                nazav(kk)=nazav(kk)+1;
                delta(kk)=dpraw16b(yy,xx)-azav(kk);
                azav(kk)=azav(kk)+delta(kk)./nazav(kk);
                m2(kk)=m2(kk)+delta(kk)*(dpraw16b(yy,xx)-azav(kk));
            end
        end
    end
    azvar=m2./nazav;
    % Normalized variance -------
    % nazvar=azvar./(azav.*azav);
    % handles.nazvar = nazvar;
    handles.azav = azav;
    handles.azvar = azvar;
    guidata(hObject,handles)
    % -------------------------------------------------------
    % Plot and save average and variance data
    % -------------------------------------------------------
    % x-axis for plots (pixel index)
    pix_xax = linspace(1,azsize,azsize);
    q_xax = pix_xax' * handles.ds*2*pi;

    figure(DPfig);
    colormap(gca,colorcube);

    rectangle('Position',[xcentre 1 nx-xcentre ycentre-1],...
        'FaceColor',[1 1 1 0.9]);

    plot(pix_xax+xcentre,ycentre-1-azav*(0.9*ycentre)/max(azav),...
        'Color', Color(3), 'DisplayName','azimuthal average [a.u.]', ...
        'linewidth', 1.0);
    plot(pix_xax+xcentre,ycentre-1-azvar*(0.9*ycentre)/max(azvar),...
        'Color', Color(11), 'DisplayName','azimuthal variance [a.u.]', ...
        'linewidth', 1.0, ...
        'LineStyle', '--');
    legend('show');
    saveas(DPfig, [pname , '\', fname{fNo}(1:end-4),'_DPCenter.fig'])

    % Azimuthal Average
    if fNo == 1
        AzimAv = figure('Name','Azimuthal Average', ...
            'NumberTitle','off');
    else
        figure(AzimAv);
        hold on
    end
    plot(q_xax,azav, 'DisplayName', fname{fNo}, 'Color', ...
        Color(3+7.5/mx0*(fNo-1)));
    xlabel('$$q$$ [\AA$$^{-1}$$]');
    ylabel('Intensity');
    legend('show', 'Interpreter', 'Latex');

    % Azimuthal Variance
    if fNo == 1
        AzimVar = figure('Name','Azimuthal Variance', ...
            'NumberTitle','off');
    else
        figure(AzimVar);
        hold on
    end
    plot(q_xax,azvar, 'DisplayName', fname{fNo}, 'Color', ...
        Color(3+7.5/mx0*(fNo-1)));
    xlabel('$$q$$ [\AA$$^{-1}$$]');
    ylabel('Variance');
    legend('show', 'Interpreter', 'Latex');

    % -------------------------------------------------------
    % Set directory to save files
    % -------------------------------------------------------
    folder = pname;
    addpath(folder);
    handles.folder = folder;
    guidata(hObject,handles)

    azav_name = sprintf('%s/%s_azav.txt',handles.folder,handles.DPfname);
    save (azav_name,'azav','-ASCII');

    azvar_name = sprintf('%s/%s_azvar.txt',handles.folder,handles.DPfname);
    save (azvar_name,'azvar','-ASCII');

    guidata(hObject,handles)
    % -------------------------------------------------------
    % compute average of all input images
    if fNo == 1
        sum_azav_name = ['_', fname{fNo}(1:5)];
        sum_azav = azav;
    else
        sum_azav_name = [sum_azav_name, '_', fname{fNo}(1:5)];
        lenDif = size(azav,1) - size(sum_azav,1);
        if lenDif > 0
            sum_azav = [sum_azav; zeros(lenDif,1)] + azav;
        elseif lenDif < 0
            sum_azav = [azav; zeros(abs(lenDif),1)] + sum_azav;
        else
            sum_azav = sum_azav + azav;
        end
    end

end
sum_azav = sum_azav./fileNumber;
sum_azav_Name = sprintf(...
    ['%s_sum_azav', sum_azav_name,'.txt'],handles.folder);
save(sum_azav_Name,'sum_azav','-ASCII');

AzimAv_Name = sprintf(...
    ['%s_AzimAv', sum_azav_name,'.fig'],handles.folder);
AzimVar_Name = sprintf(...
    ['%s_AzimVar', sum_azav_name,'.fig'],handles.folder);

saveas(AzimAv, AzimAv_Name)
saveas(AzimVar, AzimVar_Name)



%% ------------------------------------------------------------------------
% ------------------ OPEN 1D DIFFRACTION PROFILE --------------------------
% -------------------------------------------------------------------------
function Button_OpenData_Callback(hObject, eventdata, handles)

pathname = uigetdir('Test_data',...
    'Choose directory');
if pathname == 0
    % User clicked the Cancel button
    return;
end

filepath = dirrec(pathname, '_sum_azav*');  % cell array of all full directory paths including filename
if size(filepath) == [1, 0]
    disp('No files beginning with "_sum_azav" have been found.');
    return;
end
fn = filepath;
foldername = filepath;
for i = 1:size(filepath,2)
    j = 0;
    while filepath{1,i}(end-j) ~= '\'
        j = j+1;
    end
    k = j+1;
    while filepath{1,i}(end-k) ~= '\'
        k = k+1;
    end
    fn{1,i} = filepath{1,i}(end-k+1:end); % last foldername and filename e.g. '20190307_04_GeTe_105C006h_ps12a7\_sum_azav_* .txt' (for the list being prompted)
    foldername{1,i} = filepath{1,i}(end-k+1:end-j-1); % filename  e.g. '20190307_04_GeTe_105C006h_ps12a7'
    filename{1,i} = filepath{1,i}(end-j+1:end);  % filename '_sum_azav_* .txt'
end

[indx, ~] = listdlg('PromptString',{'Select the files to be processed.'},...
    'ListString',fn,'ListSize',[500 300]);

filepath   = filepath(1,indx(:));
foldername = foldername(1,indx(:));
filename   = filename(1,indx(:));
fn         = fn(1,indx(:));

if iscell(filename) ~= 1
    if filename == 0
        % User clicked the Cancel button
        return;
    end
    filename = char2cell(filename);
end

fileNumber         = size(filename,2);
handles.fileNumber = fileNumber;
handles.datfname   = foldername;
handles.datpath    = pathname;

dat = cell(1,fileNumber);
addpath(pathname);
rehash toolboxcache;


% Filter TagName from foldername
TagParts = strings(fileNumber,16);
for fNo = 1:fileNumber
    c = 1;  % count of delimiter '_' in given string
    for s = 1:length(foldername{fNo})
        if foldername{fNo}(s) ~= '_'
            TagParts(fNo,c) = TagParts(fNo,c) + foldername{fNo}(s);
        else
            c = c + 1;
        end
    end
end

TagName          = strings(fileNumber,1);
TagNameLegend    = strings(fileNumber,1);
TagNameIdentical = strings(fileNumber,1);
TagTitle = string;
identical = true;
for c = 1:size(TagParts,2)
    TagNameIdentical(:,1) = TagNameIdentical(:,1) + ' ' + TagParts(:,c);
    if ~all(TagParts(:,c) == TagParts(1,c))
        identical = false;
        for fNo = 1:fileNumber
            if TagParts(fNo,c) ~= ""
                stringg = convertStringsToChars(TagParts(fNo,c));

                % sort out probe series ID (e.g. ps17b3)
                if stringg(1:2) ~= "ps"
                    len1 = length(convertStringsToChars(TagName(fNo,1)));
                    len2 = length(convertStringsToChars(TagNameLegend(fNo,1)));

                    sk1 = strfind(stringg,'C');    % C
                    sk2 = strfind(stringg,'s');    % s
                    sk3 = strfind(stringg,'min');  % min
                    sk4 = strfind(stringg,'h');    % h
                    sk5 = strfind(stringg,'W');    % W
                    if stringg == "asd"  % convert "asd" to "as-dep."
                        TagName(fNo,1)       = TagName(fNo,1) + ' as-dep.';
                        TagNameLegend(fNo,1) = TagNameLegend(fNo,1) + ' as-dep.';
                    elseif isempty(sk1) == false  % convert "#C" to "#°C"
                        [num, status] = str2num(stringg(1:sk1-1));
                        if status == 1
                            TagName(fNo,1)       = TagName(fNo,1) + ' ' + num + ' °C ' + stringg(sk1+1:end);
                            TagNameLegend(fNo,1) = TagNameLegend(fNo,1) + ' ' + num + '$\,^{\circ}$C ' + stringg(sk1+1:end);
                        else
                            TagName(fNo,1)       = TagName(fNo,1) + ' ' + TagParts(fNo,c);
                            TagNameLegend(fNo,1) = TagNameLegend(fNo,1) + ' ' + TagParts(fNo,c);
                        end
                    else  % numbers 02, 03, 04, ...
                        TagName(fNo,1) = TagName(fNo,1) + ' ' + TagParts(fNo,c);
                        TagNameLegend(fNo,1) = TagNameLegend(fNo,1) + ' ' + TagParts(fNo,c);
                    end

                    % convert "s" to "$$\;s$$" (not in "as-dep.")
                    if isempty(sk2) == false
                        if length(stringg) == sk2 || stringg(sk2+1) ~= 'd'

                            x1 = " " + TagName{fNo,1}(1:len1);
                            y1 = strrep(TagName{fNo,1}(len1+1:end), 's', ' s');
                            TagName(fNo,1) = x1 + y1;

                            x2 = " " + TagNameLegend{fNo,1}(1:len2);
                            y2 = strrep(TagNameLegend{fNo,1}(len2+1:end), 's', '$\,$s');
                            TagNameLegend(fNo,1) = x2 + y2;

                        end
                    end

                    % convert "min" to "$$\;min$$"
                    if isempty(sk3) == false

                        x1 = " " + TagName{fNo,1}(1:len1);
                        y1 = strrep(TagName{fNo,1}(len1+1:end), 'min', ' min');
                        TagName(fNo,1) = x1 + y1;

                        x2 = " " + TagNameLegend{fNo,1}(1:len2);
                        y2 = strrep(TagNameLegend{fNo,1}(len2+1:end), 'min', '$\,$min');
                        TagNameLegend(fNo,1) = x2 + y2;

                    end

                    % convert "h" to "$$\;h$$"
                    if isempty(sk4) == false

                        x1 = " " + TagName{fNo,1}(1:len1);
                        y1 = strrep(TagName{fNo,1}(len1+1:end), 'h', ' h');
                        TagName(fNo,1) = x1 + y1;

                        x2 = " " + TagNameLegend{fNo,1}(1:len2);
                        y2 = strrep(TagNameLegend{fNo,1}(len2+1:end), 'h', '$\,$h');
                        TagNameLegend(fNo,1) = x2 + y2;
                    end

                    % convert "W" to "$$\;W$$"
                    if isempty(sk5) == false

                        x1 = " " + stringg(1:len1);
                        y1 = strrep(stringg(len1+1:end), 'W', ' W');
                        TagName(fNo,1) = x1 + y1;

                        x2 = " " + stringg(1:len2);
                        y2 = strrep(stringg(len2+1:end), 'W', '$\,$W');
                        TagNameLegend(fNo,1) = x2 + y2;
                    end
                end
            end
        end
    else
        TagTitle = TagTitle + ' ' + TagParts(1,c);
    end
end

if identical == true
    for fNo = 1:fileNumber
        TagName(fNo,1) = num2str(fNo); %+TagNameIdentical(fNo,1);
    end
end

% remove double, triple or quadruple spaces
TagName  = strrep(TagName, '    ', ' ');
TagName  = strrep(TagName, '   ', ' ');
TagName  = strrep(TagName, '  ', ' ');
TagNameLegend  = strrep(TagNameLegend, '    ', ' ');
TagNameLegend  = strrep(TagNameLegend, '   ', ' ');
TagNameLegend  = strrep(TagNameLegend, '  ', ' ');

% remove leading and trailing spaces
TagName  = strtrim(TagName);
TagName  = cellstr(TagName');
TagNameLegend  = strtrim(TagNameLegend);
TagNameLegend  = cellstr(TagNameLegend');

TagTitle = string(TagTitle);
TagTitle = strtrim(TagTitle);


% sans serif textstyle
for fNo = 1:fileNumber
    TagNameLegend{fNo} = ['\textsf{' TagNameLegend{fNo} '}'];
end

assignin('base', 'TagName', TagName);
assignin('base', 'TagNameLegend', TagNameLegend);
assignin('base', 'TagTitle', TagTitle);


% mx needed for the equidistant Colors in the plots.
if fileNumber == 1
    mx = 1;
else
    mx = fileNumber - 1;
end
handles.mx = mx;

cla(handles.axes1)
axes(handles.axes1);
for fNo = 1:fileNumber
    dat{fNo} = load(filepath{1,fNo},'-ascii');

    [pathstr,name,ext] = fileparts(filename{1,fNo}); %#ok<ASGLU>
    guidata(hObject,handles)

    % Plot data
    plot(dat{fNo}, 'DisplayName', foldername{fNo}(:),...
        'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
    legend('Interpreter', 'none');
    xlabel('Pixel Values');
    ylabel('Pixel Intensity');

    % Find maximal plot range for all curves
    if fNo == 1
        points = length(dat{fNo});
    else
        points = min(points, length(dat{fNo}));
    end

    %index number
    num = linspace(1,points,points);
    index = num.'; %transpose
end

handles.dat = dat;
handles.index = index;
set(handles.d1, 'String', 1);
set(handles.d2, 'String', points);
guidata(hObject,handles)
Button_Plot_Callback(handles.Button_Plot, eventdata, handles);

% --------------------- BEGINNING OF DATA RANGE ---------------------------
function d1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function d1_Callback(hObject, eventdata, handles)
d1 = str2double(get(hObject,'String'));
if isnan(d1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
set(hObject, 'String', d1);
guidata(hObject,handles)

% ----------------------- END OF DATA RANGE -------------------------------
function d2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function d2_Callback(hObject, eventdata, handles)
d2 = str2double(get(hObject,'String'));
if isnan(d2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if d2 > length(handles.index)
    max_index = max(handles.index);
    set(hObject, 'String', max_index);
    errordlg('Input exceeds data range','Error');
end
set(hObject, 'String', d2);
guidata(hObject,handles)

% ---------------- PLOT BEGINNING AND END OF DATA RANGE -------------------
function Button_Plot_Callback(hObject, eventdata, handles)
cla(handles.axes2);
cla(handles.axes3);

fileNumber = handles.fileNumber;
handles.I = handles.dat;
handles.Firsttime = 1;

d1 = str2double(get(handles.d1, 'String'));
d2 = str2double(get(handles.d2, 'String'));

% beginning of data range
d1b = uint16((d2-d1)/4 + d1);
% end of data range
d2a = uint16(d2 - (d2-d1)/4);

% corresponding x axes (pixel values)
x1 = handles.index(d1:d1b);
x2 = handles.index(d2a:d2);
mx = handles.mx;

% plot beginning
axes(handles.axes2);
for fNo = 1:fileNumber

    % new range of data points to be analysed
    beginning = handles.dat{fNo}(d1:d1b);

    % plot new range of data points
    legend off
    plot(x1,beginning, 'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
end
xlabel('Pixel Values');
ylabel('Pixel Intensity');
legend('Beginning of data range');
legend boxoff;

% plot ending
axes(handles.axes3);
for fNo = 1:fileNumber

    % new range of data points to be analysed
    ending = handles.dat{fNo}(d2a:d2);

    % plot new range of data points
    legend off
    plot(x2,ending, 'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
end
xlabel('Pixel Values');
ylabel('Pixel Intensity');
legend('End of data range');
legend boxoff;

% redefine selected data range
handles.x = handles.index(d1:d2);
for fNo = 1:fileNumber
    handles.I{fNo} = handles.dat{fNo}(d1:d2);
end

guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ------------------------ PLOT CALIBRATED I(Q) ---------------------------
% -------------------------------------------------------------------------
function Button_Iq_Callback(hObject, eventdata, handles)

% convert pixel values to q
index      = handles.x;
fileNumber = handles.fileNumber;
TagName       = evalin('base', 'TagName');
TagNameLegend = evalin('base', 'TagNameLegend');
mx = handles.mx;

q = index*handles.ds*2*pi;


% % input q data for GeTe_Lamellae
% d1  = str2double(get(handles.d1, 'String'));
% d2  = str2double(get(handles.d2, 'String'));
%
% filep = 'C:\Users\aqchr\Documents\Masterarbeit\01_Daten_von_Geraeten\03_Data_TEM\02_SAED\GeTe_Lamella\';
%
% q_asd       = load([filep, '_sum_azav_24pat_asd_x.txt'], '-ascii');
% q_asdtecnai = load([filep, '_sum_azav_24pat_asd_tecnai_x.txt'], '-ascii');
%
% q   = q_asdtecnai(d1:d2)/10*(2*pi);


handles.q = q;

I = handles.I;
%-------------------------- Normalization ---------------------------------
% Find maximum of intensity I for each curve
Imax = cell(fileNumber,1);
for fNo = 1:fileNumber
    Imax{fNo} = max(I{fNo});
end
Imax = cell2mat(Imax);

% Scale all curves
Imaxall = max(Imax(:));
for fNo = 1:fileNumber
    I{fNo} = I{fNo}/Imax(fNo)*Imaxall;
end
handles.I = I;
% -------------------------------------------------------------------------

% Plot I(q) vs q
if handles.Firsttime == 1
    cla(handles.axes4, 'reset');
    axes(handles.axes4);
    for fNo = 1:fileNumber
        plot(q,I{fNo}, 'DisplayName', TagNameLegend{fNo}(:)',...
            'Color', Color(3+7.5/mx*(fNo-1)));
        hold on
        legend('Interpreter', 'Latex');
        xlabel('$$q$$ [\AA$$^{-1}$$]');
        ylabel('$$I(q)$$');
    end
    handles.Firsttime = 0;
end

% convert q to s
s = q/2/pi;
handles.s = s;

% s^2
s2 = s.^2;
handles.s2 = s2;

% Length of s/q
L = uint16(length(q));
handles.L = L;

guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ------------------------ ELEMENT HANDLES --------------------------------
% -------------------------------------------------------------------------

% --------------------- ELEMENT SELECT POPUP ------------------------------
function Element1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Element1_Callback(hObject, eventdata, handles)

% get selected element index / atomic number
val = get(hObject,'Value');
elem1 = val - 1;
handles.elem1 = elem1;
% get content string for output
contents = cellstr(get(hObject,'String'));
EName1 = contents{val};
handles.EName1 = EName1;

guidata(hObject,handles)

function Element2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Element2_Callback(hObject, eventdata, handles)
% get selected element index / atomic number
val = get(hObject,'Value');
elem2 = val - 1;
handles.elem2 = elem2;
% get content string for output
contents = cellstr(get(hObject,'String'));
EName2 = contents{val};
handles.EName2 = EName2;

guidata(hObject,handles)

function Element3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Element3_Callback(hObject, eventdata, handles)
% get selected element index / atomic number
val = get(hObject,'Value');
elem3 = val - 1;
handles.elem3 = elem3;
% get content string for output
contents = cellstr(get(hObject,'String'));
EName3 = contents{val};
handles.EName3 = EName3;

guidata(hObject,handles)

function Element4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Element4_Callback(hObject, eventdata, handles)
% get selected element index / atomic number
val = get(hObject,'Value');
elem4 = val - 1;
handles.elem4 = elem4;
% get content string for output
contents = cellstr(get(hObject,'String'));
EName4 = contents{val};
handles.EName4 = EName4;

guidata(hObject,handles)

function Element5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Element5_Callback(hObject, eventdata, handles)
% get selected element index / atomic number
val = get(hObject,'Value');
elem5 = val - 1;
handles.elem5 = elem5;
% get content string for output
contents = cellstr(get(hObject,'String'));
EName5 = contents{val};
handles.EName5 = EName5;

guidata(hObject,handles)

% ------------------------ ELEMENT FRACTIONS ------------------------------
function edit_e1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_e1_Callback(hObject, eventdata, handles)
% --- Composition of element 1
e1 = str2double(get(hObject,'String'));
if isnan(e1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.e1 = e1;

guidata(hObject,handles)

function edit_e2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_e2_Callback(hObject, eventdata, handles)
% --- Composition of element 2
e2 = str2double(get(hObject,'String'));
if isnan(e2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.e2 = e2;

guidata(hObject,handles)

function edit_e3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_e3_Callback(hObject, eventdata, handles)
% --- Composition of element 3
e3 = str2double(get(hObject,'String'));
if isnan(e3)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.e3 = e3;

guidata(hObject,handles)

function edit_e4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_e4_Callback(hObject, eventdata, handles)
% --- Composition of element 4
e4 = str2double(get(hObject,'String'));
if isnan(e4)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.e4 = e4;

guidata(hObject,handles)

function edit_e5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_e5_Callback(hObject, eventdata, handles)
% --- Composition of element 5
e5 = str2double(get(hObject,'String'));
if isnan(e5)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.e5 = e5;

guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ------------------------ FITTING HANDLES --------------------------------
% -------------------------------------------------------------------------

% --------------------------- FIT RANGE -----------------------------------
function slider_fitrange_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_fitrange_Callback(hObject, eventdata, handles)
% Fit only over tail end (last (Perc)% of q range)
Perc = round(get(hObject, 'Value'));
handles.Perc = Perc;
set(handles.text_fitrange,'String',['Last ',num2str(Perc),'% of curve']);

guidata(hObject,handles)


% ------------------------ CURVE SELECTION --------------------------------
function popSelectCurve_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', 'No Data');
assignin('base','SelectedFitCurve', get(hObject, 'Value'));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popSelectCurve_Callback(hObject, eventdata, handles)
SFC = get(hObject, 'Value');
assignin('base','SelectedFitCurve', SFC);
ResNAll     = evalin('base','ResNAll');
Resq_fitAll = evalin('base','Resq_fitAll');
ResSSAll    = evalin('base','ResSSAll');
ResDamp     = evalin('base','ResDamp');
CDegreeM    = evalin('base','CDegreeM');


set(handles.text_CDegree, 'String', CDegreeM(SFC));
set(handles.text_SS, 'String', ResSSAll{SFC});
set(handles.text_q_fit, 'String', Resq_fitAll{SFC});
set(handles.edit_q_fit, 'String', handles.q_fixAll{SFC});
set(handles.text_N, 'String', ResNAll{SFC});
set(handles.edit_N, 'String', handles.NAll{SFC});
set(handles.text_damping, 'String', ResDamp{SFC});

% ------------------- KIRKLAND/LOBATO SELECTION ---------------------------
function popup_param_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popup_param_Callback(hObject, eventdata, handles)
% --- choose which Parameterisation factors to use in fitting
param_val = get(hObject,'Value');
switch param_val
    case 3
        paramL = load('Lobato_2014.txt','-ascii');
        handles.paramL = paramL;
    otherwise
        paramK = load('Kirkland_2010.txt','-ascii'); % default
        handles.paramK = paramK;
end

handles.param_val = param_val;
guidata(hObject,handles)

% ------------------------------ N ----------------------------------------
function edit_DeltaN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_DeltaN_Callback(hObject, eventdata, handles)
DeltaN = str2double(get(hObject,'String'));
if isnan(DeltaN)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

function edit_N_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_N_Callback(hObject, eventdata, handles)
N = str2double(get(hObject,'String'));
if isnan(N)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.N = N;
SFC = evalin('base', 'SelectedFitCurve');
handles.NAll{SFC} = N;
% print value of N in editable text  (edit_N)
% set(handles.edit_N, 'String', handles.NAll{SFC});

guidata(hObject,handles)

function edit_dN_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '0.3')
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_dN_Callback(hObject, eventdata, handles)
dN = str2double(get(hObject,'String'));
if isnan(dN)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.dN = dN;

guidata(hObject,handles)

function Nplus_Callback(hObject, eventdata, handles)
handles.N = handles.N + handles.dN;
SFC = evalin('base', 'SelectedFitCurve');
handles.NAll{SFC} = handles.NAll{SFC} + handles.dN;
% print new value of N in editable text (edit_N)
set(handles.edit_N, 'String', handles.NAll{SFC});

guidata(hObject,handles)

function Nminus_Callback(hObject, eventdata, handles)
handles.N = handles.N - handles.dN;
SFC = evalin('base', 'SelectedFitCurve');
handles.NAll{SFC} = handles.NAll{SFC} - handles.dN;
% print new value of N in editable text (edit_N)
set(handles.edit_N, 'String', handles.NAll{SFC});

guidata(hObject,handles)

% ------------------------------ q ----------------------------------------
function edit_q_fit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_q_fit_Callback(hObject, eventdata, handles)
% --- q_fix = value close to which user wants fitting to be done
q_fix = str2double(get(hObject,'String'));
if isnan(q_fix)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.q_fix = q_fix;
SFC = evalin('base', 'SelectedFitCurve');
handles.q_fixAll{SFC} = q_fix;
% print desired value of q in editable text (edit_q_fit)
% actual value (data point) will update when 'Fit Data' button is pushed
% set(handles.edit_q_fit, 'String', handles.q_fix);

guidata(hObject,handles)

function edit_dq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_dq_Callback(hObject, eventdata, handles)
% --- value to change edit_q_fit by
dq = str2double(get(hObject,'String'));
if isnan(dq)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.dq = dq;

guidata(hObject,handles)

function qplus_Callback(hObject, eventdata, handles)
% --- increase edit_q_fit by given edit_dq value
handles.q_fix = handles.q_fix + handles.dq;
SFC = evalin('base', 'SelectedFitCurve');
handles.q_fixAll{SFC} = handles.q_fixAll{SFC} + handles.dq;
% print desired value of q in static text
% actual value (data point) will update when 'Fit Data' button is pushed
set(handles.edit_q_fit, 'String', handles.q_fixAll{SFC});

guidata(hObject,handles)

function qminus_Callback(hObject, eventdata, handles)
% --- decrease edit_q_fit by given edit_dq value
handles.q_fix = handles.q_fix - handles.dq;
SFC = evalin('base', 'SelectedFitCurve');
handles.q_fixAll{SFC} = handles.q_fixAll{SFC} - handles.dq;
% print desired value of q in static text
% actual value (data point) will update when 'Fit Data' button is pushed
set(handles.edit_q_fit, 'String', handles.q_fixAll{SFC});

guidata(hObject,handles)

% ------------------------- EDIT_DAMPING ----------------------------------
function edit_damping_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_damping_Callback(hObject, eventdata, handles)
% --- set edit_damping factor (default 0.3)
damping = str2double(get(hObject,'String'));
if isnan(damping)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.edit_damping = damping;
% print value of edit_damping factor in editable text
% set(handles.edit_damping, 'String', damping);

guidata(hObject,handles)

% -------------------------- PLOT RANGE -----------------------------------
function rmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rmax_Callback(hObject, eventdata, handles)
% --- rmax sets the range of r (x-axis) to plot G(r)
rmax = str2double(get(hObject,'String'));
if isnan(rmax)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.rmax = rmax;

guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ------------------------ AUTOFIT ALL ------------------------------------
% -------------------------------------------------------------------------
function Button_Autofit_all_Callback(hObject, eventdata, handles)
%% Compute gq = <f^2(q)>
s2 = handles.s2;
elem1 = handles.elem1;
elem2 = handles.elem2;
elem3 = handles.elem3;
elem4 = handles.elem4;
elem5 = handles.elem5;
e1 = handles.e1;
e2 = handles.e2;
e3 = handles.e3;
e4 = handles.e4;
e5 = handles.e5;
% compute atomic ratio (composition)
e_tot = e1+e2+e3+e4+e5;
e_r1 = e1/e_tot;
e_r2 = e2/e_tot;
e_r3 = e3/e_tot;
e_r4 = e4/e_tot;
e_r5 = e5/e_tot;
handles.e_tot = e_tot;
handles.e_r1 = e_r1;
handles.e_r2 = e_r2;
handles.e_r3 = e_r3;
handles.e_r4 = e_r4;
handles.e_r5 = e_r5;
guidata(hObject,handles)

if handles.param_val == 3 % Lobato

    paramL = handles.paramL;
    paramL_1 = paramL(elem1,:);
    paramL_2 = paramL(elem2,:);
    paramL_3 = paramL(elem3,:);
    paramL_4 = paramL(elem4,:);
    paramL_5 = paramL(elem5,:);
    paramL_elem = [paramL_1;paramL_2;paramL_3;paramL_4;paramL_5];
    paramL_table = array2table(paramL_elem,...
        'VariableNames',{'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5'},...
        'RowNames',{'elem1','elem2','elem3','elem4','elem5'});
    handles.paramL_table = paramL_table;

    A1_1 = paramL_table{'elem1','A1'};
    A2_1 = paramL_table{'elem1','A2'};
    A3_1 = paramL_table{'elem1','A3'};
    A4_1 = paramL_table{'elem1','A4'};
    A5_1 = paramL_table{'elem1','A5'};
    B1_1 = paramL_table{'elem1','B1'};
    B2_1 = paramL_table{'elem1','B2'};
    B3_1 = paramL_table{'elem1','B3'};
    B4_1 = paramL_table{'elem1','B4'};
    B5_1 = paramL_table{'elem1','B5'};

    f1 = ((s2*B1_1+1).^2).\(A1_1*(s2*B1_1+2))+((s2*B2_1+1).^2).\(A2_1*(s2*B2_1+2))+((s2*B3_1+1).^2).\(A3_1*(s2*B3_1+2))+((s2*B4_1+1).^2).\(A4_1*(s2*B4_1+2))+((s2*B5_1+1).^2).\(A5_1*(s2*B5_1+2));

    A1_2 = paramL_table{'elem2','A1'};
    A2_2 = paramL_table{'elem2','A2'};
    A3_2 = paramL_table{'elem2','A3'};
    A4_2 = paramL_table{'elem2','A4'};
    A5_2 = paramL_table{'elem2','A5'};
    B1_2 = paramL_table{'elem2','B1'};
    B2_2 = paramL_table{'elem2','B2'};
    B3_2 = paramL_table{'elem2','B3'};
    B4_2 = paramL_table{'elem2','B4'};
    B5_2 = paramL_table{'elem2','B5'};

    f2 = ((s2*B1_2+1).^2).\(A1_2*(s2*B1_2+2))+((s2*B2_2+1).^2).\(A2_2*(s2*B2_2+2))+((s2*B3_2+1).^2).\(A3_2*(s2*B3_2+2))+((s2*B4_2+1).^2).\(A4_2*(s2*B4_2+2))+((s2*B5_2+1).^2).\(A5_2*(s2*B5_2+2));

    A1_3 = paramL_table{'elem3','A1'};
    A2_3 = paramL_table{'elem3','A2'};
    A3_3 = paramL_table{'elem3','A3'};
    A4_3 = paramL_table{'elem3','A4'};
    A5_3 = paramL_table{'elem3','A5'};
    B1_3 = paramL_table{'elem3','B1'};
    B2_3 = paramL_table{'elem3','B2'};
    B3_3 = paramL_table{'elem3','B3'};
    B4_3 = paramL_table{'elem3','B4'};
    B5_3 = paramL_table{'elem3','B5'};

    f3 = ((s2*B1_3+1).^2).\(A1_3*(s2*B1_3+2))+((s2*B2_3+1).^2).\(A2_3*(s2*B2_3+2))+((s2*B3_3+1).^2).\(A3_3*(s2*B3_3+2))+((s2*B4_3+1).^2).\(A4_3*(s2*B4_3+2))+((s2*B5_3+1).^2).\(A5_3*(s2*B5_3+2));

    A1_4 = paramL_table{'elem4','A1'};
    A2_4 = paramL_table{'elem4','A2'};
    A3_4 = paramL_table{'elem4','A3'};
    A4_4 = paramL_table{'elem4','A4'};
    A5_4 = paramL_table{'elem4','A5'};
    B1_4 = paramL_table{'elem4','B1'};
    B2_4 = paramL_table{'elem4','B2'};
    B3_4 = paramL_table{'elem4','B3'};
    B4_4 = paramL_table{'elem4','B4'};
    B5_4 = paramL_table{'elem4','B5'};

    f4 = ((s2*B1_4+1).^2).\(A1_4*(s2*B1_4+2))+((s2*B2_4+1).^2).\(A2_4*(s2*B2_4+2))+((s2*B3_4+1).^2).\(A3_4*(s2*B3_4+2))+((s2*B4_4+1).^2).\(A4_4*(s2*B4_4+2))+((s2*B5_4+1).^2).\(A5_4*(s2*B5_4+2));

    A1_5 = paramL_table{'elem5','A1'};
    A2_5 = paramL_table{'elem5','A2'};
    A3_5 = paramL_table{'elem5','A3'};
    A4_5 = paramL_table{'elem5','A4'};
    A5_5 = paramL_table{'elem5','A5'};
    B1_5 = paramL_table{'elem5','B1'};
    B2_5 = paramL_table{'elem5','B2'};
    B3_5 = paramL_table{'elem5','B3'};
    B4_5 = paramL_table{'elem5','B4'};
    B5_5 = paramL_table{'elem5','B5'};

    f5 = ((s2*B1_5+1).^2).\(A1_5*(s2*B1_5+2))+((s2*B2_5+1).^2).\(A2_5*(s2*B2_5+2))+((s2*B3_5+1).^2).\(A3_5*(s2*B3_5+2))+((s2*B4_5+1).^2).\(A4_5*(s2*B4_5+2))+((s2*B5_5+1).^2).\(A5_5*(s2*B5_5+2));

else % Kirkland

    paramK = handles.paramK;
    paramK_1 = paramK(elem1,:);
    paramK_2 = paramK(elem2,:);
    paramK_3 = paramK(elem3,:);
    paramK_4 = paramK(elem4,:);
    paramK_5 = paramK(elem5,:);
    paramK_elem = [paramK_1;paramK_2;paramK_3;paramK_4;paramK_5];
    paramK_table = array2table(paramK_elem,...
        'VariableNames',{'a1','b1','a2','b2','a3','b3','c1','d1','c2','d2','c3','d3'},...
        'RowNames',{'elem1','elem2','elem3','elem4','elem5'});
    handles.paramK_table = paramK_table;

    a1_1 = paramK_table{'elem1','a1'};
    a2_1 = paramK_table{'elem1','a2'};
    a3_1 = paramK_table{'elem1','a3'};
    b1_1 = paramK_table{'elem1','b1'};
    b2_1 = paramK_table{'elem1','b2'};
    b3_1 = paramK_table{'elem1','b3'};
    c1_1 = paramK_table{'elem1','c1'};
    c2_1 = paramK_table{'elem1','c2'};
    c3_1 = paramK_table{'elem1','c3'};
    d1_1 = paramK_table{'elem1','d1'};
    d2_1 = paramK_table{'elem1','d2'};
    d3_1 = paramK_table{'elem1','d3'};

    f1 = ((s2+b1_1).\a1_1)+((s2+b2_1).\a2_1)+((s2+b3_1).\a3_1)+(exp(-s2.*d1_1).*c1_1)+(exp(-s2.*d2_1).*c2_1)+(exp(-s2.*d3_1).*c3_1);

    a1_2 = paramK_table{'elem2','a1'};
    a2_2 = paramK_table{'elem2','a2'};
    a3_2 = paramK_table{'elem2','a3'};
    b1_2 = paramK_table{'elem2','b1'};
    b2_2 = paramK_table{'elem2','b2'};
    b3_2 = paramK_table{'elem2','b3'};
    c1_2 = paramK_table{'elem2','c1'};
    c2_2 = paramK_table{'elem2','c2'};
    c3_2 = paramK_table{'elem2','c3'};
    d1_2 = paramK_table{'elem2','d1'};
    d2_2 = paramK_table{'elem2','d2'};
    d3_2 = paramK_table{'elem2','d3'};

    f2 = ((s2+b1_2).\a1_2)+((s2+b2_2).\a2_2)+((s2+b3_2).\a3_2)+(exp(-s2.*d1_2).*c1_2)+(exp(-s2.*d2_2).*c2_2)+(exp(-s2.*d3_2).*c3_2);

    a1_3 = paramK_table{'elem3','a1'};
    a2_3 = paramK_table{'elem3','a2'};
    a3_3 = paramK_table{'elem3','a3'};
    b1_3 = paramK_table{'elem3','b1'};
    b2_3 = paramK_table{'elem3','b2'};
    b3_3 = paramK_table{'elem3','b3'};
    c1_3 = paramK_table{'elem3','c1'};
    c2_3 = paramK_table{'elem3','c2'};
    c3_3 = paramK_table{'elem3','c3'};
    d1_3 = paramK_table{'elem3','d1'};
    d2_3 = paramK_table{'elem3','d2'};
    d3_3 = paramK_table{'elem3','d3'};

    f3 = ((s2+b1_3).\a1_3)+((s2+b2_3).\a2_3)+((s2+b3_3).\a3_3)+(exp(-s2.*d1_3).*c1_3)+(exp(-s2.*d2_3).*c2_3)+(exp(-s2.*d3_3).*c3_3);

    a1_4 = paramK_table{'elem4','a1'};
    a2_4 = paramK_table{'elem4','a2'};
    a3_4 = paramK_table{'elem4','a3'};
    b1_4 = paramK_table{'elem4','b1'};
    b2_4 = paramK_table{'elem4','b2'};
    b3_4 = paramK_table{'elem4','b3'};
    c1_4 = paramK_table{'elem4','c1'};
    c2_4 = paramK_table{'elem4','c2'};
    c3_4 = paramK_table{'elem4','c3'};
    d1_4 = paramK_table{'elem4','d1'};
    d2_4 = paramK_table{'elem4','d2'};
    d3_4 = paramK_table{'elem4','d3'};

    f4 = ((s2+b1_4).\a1_4)+((s2+b2_4).\a2_4)+((s2+b3_4).\a3_4)+(exp(-s2.*d1_4).*c1_4)+(exp(-s2.*d2_4).*c2_4)+(exp(-s2.*d3_4).*c3_4);

    a1_5 = paramK_table{'elem5','a1'};
    a2_5 = paramK_table{'elem5','a2'};
    a3_5 = paramK_table{'elem5','a3'};
    b1_5 = paramK_table{'elem5','b1'};
    b2_5 = paramK_table{'elem5','b2'};
    b3_5 = paramK_table{'elem5','b3'};
    c1_5 = paramK_table{'elem5','c1'};
    c2_5 = paramK_table{'elem5','c2'};
    c3_5 = paramK_table{'elem5','c3'};
    d1_5 = paramK_table{'elem5','d1'};
    d2_5 = paramK_table{'elem5','d2'};
    d3_5 = paramK_table{'elem5','d3'};

    f5 = ((s2+b1_5).\a1_5)+((s2+b2_5).\a2_5)+((s2+b3_5).\a3_5)+(exp(-s2.*d1_5).*c1_5)+(exp(-s2.*d2_5).*c2_5)+(exp(-s2.*d3_5).*c3_5);
end

fq = (f1.*e_r1) + (f2.*e_r2) + (f3.*e_r3) + (f4.*e_r4) + (f5.*e_r5);
fq_sq = fq.^2;
handles.fq_sq = fq_sq;

%% Compute gq = <f^2(q)>
gq = (f1.^2*e_r1) + (f2.^2*e_r2) + (f3.^2*e_r3) + (f4.^2*e_r4) + (f5.^2*e_r5);
handles.gq = gq;
guidata(hObject,handles)

%% Auto Fit atomic scattering curve = N*gq+C
Iq = handles.I;
gq = handles.gq;
q = handles.q;
L = handles.L;
fileNumber = handles.fileNumber;
mx = handles.mx;
TagName       = evalin('base', 'TagName');
TagNameLegend = evalin('base', 'TagNameLegend');
TagTitle      = evalin('base', 'TagTitle');
SFC           = evalin('base', 'SelectedFitCurve');


% Correct high-frequency noise of Iq for q>smind (mean interpolation)

Iqq       = Iq;
IqqSmooth = Iq;
smind = 1000;  % index where smoothness is increased
for fNo = 1:fileNumber
    Iqq{fNo}                  = Iq{fNo};
    IqqSmooth{fNo}(1:smind)   = smooth(Iqq{fNo}(1:smind), 3);
    IqqSmooth{fNo}(smind:end) = smooth(Iqq{fNo}(smind:end), 35);
    Iq{fNo}                   = IqqSmooth{fNo};
end
handles.I = Iq;

% Selection of weights based on Fit range
% Fit only over tail end (last (Perc)% of q range)
Perc = handles.Perc;
AFrange = 0.01*(100 - Perc)*handles.L;
wi = ones(L,1);
wi(1:AFrange) = 0;

% qmax: point where fitted curve crosses Iq curve
[qmax,qpos] = max(q); %value and array position of qmax
fqfit = gq(qpos);

handles.AutofitAll = handles.I;
handles.NAll = handles.I;
for fNo = 1:fileNumber

    iqfit = sum(Iq{fNo}(qpos));

    a1 = sum(wi.*gq.*Iq{fNo});
    a2 = sum(wi.*Iq{fNo}*fqfit);
    a3 = sum(wi.*gq*iqfit);
    a4 = sum(wi)*fqfit*iqfit;
    a5 = sum(wi.*gq.^2);
    a6 = 2*sum(wi.*gq*fqfit);
    a7 = sum(wi)*fqfit*fqfit;

    N = (a1-a2-a3+a4)/(a5-a6+a7);
    C = iqfit-N*fqfit;

    Autofit = N*gq+C;
    handles.AutofitAll{fNo} = Autofit;
    handles.NAll{fNo}       = N;
    guidata(hObject,handles)
end


% Plot auto fit curve
cla(handles.axes4);
axes(handles.axes4);
FitAx4 = cell(1,fileNumber);
for fNo = 1:fileNumber
    plot(handles.q, handles.I{fNo}, 'DisplayName', [TagNameLegend{fNo}(:)' '\textsf{ + Fit}'],...
        'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
    FitAx4{fNo} = plot(handles.q, handles.AutofitAll{fNo}, '--',...
        'Color', Color(3+7.5/mx*(fNo-1)), 'Handlevisibility', 'off');
end
xlabel('$$q$$ [\AA$$^{-1}$$]');
ylabel('$$I(q)$$');
legend('show', 'Interpreter', 'Latex');
assignin('base','FitAx4',FitAx4);

% Plot magnified view
cla(handles.axes5, 'reset');
axes(handles.axes5);
FitAx5 = cell(1,fileNumber);
q2 = handles.q(L/2:L);
handles.q2 = q2;
for fNo = 1:fileNumber
    handles.mag_fit = handles.AutofitAll{fNo}(L/2:L);
    handles.mag_dat = handles.I{fNo}(L/2:L);
    plot(handles.q2, handles.mag_dat, 'DisplayName', [TagNameLegend{fNo}(:)' '\textsf{ + Fit}'], ...
        'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
    FitAx5{fNo} = plot(handles.q2, handles.mag_fit, '--', ...
        'Color', Color(3+7.5/mx*(fNo-1)), 'Handlevisibility', 'off');
end
xlabel('$$q$$ [\AA$$^{-1}$$]');
ylabel('$$I(q)$$');
legend('show', 'Interpreter', 'Latex');
assignin('base','FitAx5',FitAx5);

% Update GUI display (shows result of selected curve)
ResSSAll = handles.I;
ResDamp = handles.I;
handles.q_fixAll = handles.I;
for fNo = 1:fileNumber
    % goodness of fit (sum of squared deviations)
    ResSSAll{fNo} = sprintf('%12.6e',...
        sum((handles.I{fNo}-handles.AutofitAll{fNo}).^2));
    ResDamp{fNo} = handles.edit_damping;
    handles.q_fixAll{fNo} = qmax;
end

% update values of parameters on GUI
set(handles.text_SS, 'String', ResSSAll{SFC});
set(handles.text_N, 'String', handles.NAll{SFC});
set(handles.edit_N, 'String', handles.NAll{SFC});
set(handles.text_q_fit, 'String', handles.q_fixAll{SFC});
set(handles.edit_q_fit, 'String', handles.q_fixAll{SFC});
set(handles.text_damping, 'String', handles.edit_damping);

assignin('base','ResSSAll',ResSSAll);
assignin('base','ResDamp',ResDamp);
guidata(hObject,handles)
% ----------------------------------------------------------------------
q = handles.q;
I = handles.I;
I = Iq;
s = handles.s;
s2 = handles.s2;
ds = handles.ds;
fq_sq = handles.fq_sq;
NAll = handles.NAll;
AutofitAll = handles.AutofitAll;
damping = handles.edit_damping;

% Compute phiq (uncorrected)
phiq = I;
for fNo = 1:fileNumber
    phiq{fNo} = ((I{fNo} - AutofitAll{fNo}).*s)./(NAll{fNo}*fq_sq);
end

% Correct high-frequency noise of phiq for q>qnoise (mean interpolation)
qnoise = 3.;  % beginning of smoothing
wexp   = 2;
[~,iqnoise] = min(abs(q-qnoise));

for fNo = 1:fileNumber
    iFTphi  = smooth(phiq{fNo}, 35);
    phiqHfC = phiq{fNo};
    for iq = iqnoise:length(q)
        phiqHfC(iq) = (phiq{fNo}(iq)+iFTphi(iq)*(q(iq)-qnoise)^wexp)/ ...
            (1+(q(iq)-qnoise)^wexp);
    end
    phiq{fNo} = phiqHfC;
end

assignin('base','phiq', phiq);

%% Correct phiq by subtracting back-FT of G(r<maxCorrRange)
rmax = handles.rmax;
dr   = 0.01;
r    = 0.01:dr:rmax;
numMean       = 2;  %(2) number of correction ranges over which is averaged: rcoind-numMean:rcoind+numMean
CorrSteps     = 15;  %(15)
maxCorrRange  = 1.4;  %(1.4) in Angstrom
CStep         = maxCorrRange/(CorrSteps-1);  % in Angstrom
rc            = round(((1:CorrSteps)-1)*CStep,2);  % correction ranges in Angstrom

rco           = 1.0;  %(1.0) in Angstrom, Correction Range around which is averaged
[~,rcoind]    = min(abs(rco-rc));  % index with rc(index) closest to rco
MeanRangeInd  = rcoind-numMean:rcoind+numMean;  % index range over which is averaged

assignin('base', 'dr', dr);
assignin('base', 'r', r);
assignin('base', 'CorrSteps', CorrSteps);
assignin('base', 'numMean', numMean);
assignin('base', 'maxCorrRange', maxCorrRange);
assignin('base', 'rc', rc);
assignin('base', 'rco', rco);
assignin('base', 'rcoind', rcoind);
assignin('base', 'MeanRangeInd', MeanRangeInd);

FTphi     = cell(1,fileNumber);
f_cf      = cell(CorrSteps,fileNumber);
phiqCor   = cell(CorrSteps,fileNumber);
Sq        = cell(CorrSteps,fileNumber);
phiq_damp = cell(CorrSteps,fileNumber);

f_cfM      = cell(1,fileNumber);  % correction function f_{cf}
phiqCorM   = cell(1,fileNumber);
SqM        = cell(1,fileNumber);
phiq_dampM = cell(1,fileNumber);
IqCor      = cell(1,fileNumber);
FqCor      = cell(1,fileNumber);
FqNCor     = cell(1,fileNumber);

Lorch = sin(pi*q./qmax)./(pi*q./qmax);

%damping_function = exp(-s2'.*damping);
damping_function = Lorch';

for fNo = 1:fileNumber
    FTphi{fNo} = 8 * pi * phiq{fNo}'* sin(q*r) * ds; % * 3.0;
    FTSys      = 0*FTphi{fNo};

    % iC=1 -> no Correction, iC=CorrSteps -> max Correction
    for iC = 1:CorrSteps
        CorrRange        = 1:round((((iC-1)*CStep)/dr));  % in r-axis units (dr)
        FTSys(CorrRange) = FTphi{fNo}(CorrRange);
        f_cf{iC,fNo}     = 1./(2*pi) * FTSys * sin(q*r)' * dr ;
        phiqCor{iC,fNo}  = phiq{fNo}' - f_cf{iC,fNo};
        Sq{iC,fNo}       = phiqCor{iC,fNo}./q'+1;

        phiq_damp{iC,fNo} = phiqCor{iC,fNo}.*damping_function;
    end

    % compute averages over numMean largest correction ranges
    f_cfM{fNo}      = mean(cell2mat(f_cf(MeanRangeInd, fNo)));
    phiqCorM{fNo}   = mean(cell2mat(phiqCor(MeanRangeInd, fNo)));
    SqM{fNo}        = mean(cell2mat(Sq(MeanRangeInd, fNo)));
    phiq_dampM{fNo} = mean(cell2mat(phiq_damp(MeanRangeInd, fNo)));

    % compute corrected Bragg intensity (I(q) corrected)
    IqCor{fNo} = AutofitAll{fNo}' + (phiqCorM{fNo}*NAll{fNo}.*fq_sq')./s';
    IqCor{fNo} = IqCor{fNo}';

    % compute corrected F(q) and normalized F(q) for RMC fit (notation as in RMCProfile manual)
    %FqCor{fNo}  = fq_sq.*(SqM{fNo}'-1);
    FqCor{fNo}  = gq.*(SqM{fNo}'-1);
    FqNCor{fNo} = SqM{fNo}'-1;

end

%% Compute Correction Degree (CDegree)
A = zeros(CorrSteps,fileNumber);
B = zeros(1,fileNumber);  % Normalization
for fNo = 1:fileNumber
    for iC = 1:CorrSteps
        A(iC,fNo) = trapz(q, abs(f_cf{iC,fNo}));
    end
    B(1,fNo) = trapz(q, abs(phiq{fNo}));
end

CDegree   = A./B;
CDegreeM  = mean(CDegree(MeanRangeInd,:),1);
CDegreeME = std(CDegree(MeanRangeInd,:),1);
CDegreeMErel = CDegreeME./CDegreeM;

%% Update workspace
assignin('base', 'FTphi', FTphi);

assignin('base', 'f_cf', f_cf);
assignin('base', 'phiqCor', phiqCor);
assignin('base', 'Sq', Sq);
assignin('base', 'phiq_damp', phiq_damp);

assignin('base', 'f_cfM', f_cfM);
assignin('base', 'phiqCorM', phiqCorM);
assignin('base', 'SqM', SqM);
assignin('base', 'phiq_dampM', phiq_dampM);

assignin('base', 'IqCor', IqCor);
assignin('base', 'FqCor', FqCor);
assignin('base', 'FqNCor', FqNCor);

assignin('base', 'A', A);
assignin('base', 'B', B);
assignin('base', 'CDegree', CDegree);
assignin('base', 'CDegreeM', CDegreeM);
assignin('base', 'CDegreeME', CDegreeME);

% Update GUI display for CDegree (shows result of selected curve)
set(handles.text_CDegree, 'String', CDegreeM(SFC));
guidata(hObject,handles)

% create file with meta-data (saving file, when exporting Iq, phiy, Gr)
d1 = str2double(get(handles.d1, 'String'));
d2 = str2double(get(handles.d2, 'String'));
metadata = {ds; [d1,d2]; Perc; damping; q'; r; rco; maxCorrRange; numMean; CorrSteps};
assignin('base', 'metadata', metadata);

%% Plot phiq etc.
cla(handles.axes6);
axes(handles.axes6);
FitAx6a = cell(1,fileNumber);  % phi_q
FitAx6b = cell(1,fileNumber);  % phi_q corrected
FitAx6c = cell(1,fileNumber);  % phi_q_damped
FitAx6d = cell(1,fileNumber);  % diff
for fNo = 1:fileNumber
    % Plot phiq (uncorrected)
    FitAx6a{fNo} = plot(q, phiq{fNo}, 'Linestyle', '--', ...'Color', 'black');
        'Color', Color(3+7.5/mx*(fNo-1)));
    hold on

    % Plot inverse FT of G(r<maxCorrRange) to see what is subtracted
    FitAx6d{fNo} = plot(q, f_cfM{fNo}, 'k', 'Linestyle', '-.');

    % Plot corrected phiq
    FitAx6b{fNo} = plot(q, phiqCorM{fNo}, '-', ...
        'Color', Color(3+7.5/mx*(fNo-1)));

    %     % Plot corrceted phiq with damping
    %     FitAx6c{fNo} = plot(q, phiq_damp{fNo}, '--', ...
    %         'Color', Color(3+7.5/mx*(fNo-1)));

end
xlabel('$$q$$ [\AA$$^{-1}$$]');
ylabel('$$\phi(q)$$');
legend('show',  '\textsf{$$\phi_\mathrm{u}(q)$$ (uncorrected)}', ...
    '\textsf{$$f_\mathrm{cf}(q)$$ (correction function)}', ...
    '\textsf{$$\phi(q)$$ (corrected)}', ... %'\phi(q) \cdot exp(-bq^2)', ...
    'Interpreter', 'latex', 'Location', 'Southeast');


% legend updated below (after inverse FT)
assignin('base','FitAx6a',FitAx6a);
assignin('base','FitAx6b',FitAx6b);
assignin('base','FitAx6c',FitAx6c);
assignin('base','FitAx6d',FitAx6d);

% plot reference line at phiq=zero
xlims = get(handles.axes6, 'xlim');
hold on
plot([0 xlims(2)], [0 0],'k', 'HandleVisibility','off');
hold off

% ----------------------------------------------------------------------
%%
% figure('Name', 'Corrected I(q)', ...
%     'units', 'normalized', ...
%     'outerposition', [0.1 0 0.7 1], ...
%     'NumberTitle', 'off')
%
% % Plot corrected I(q) (for RMC refinement)
% plot(q, I{1}, 'Linestyle', '--', ... 'Color', 'black');
%     'Color', Color(3+7.5/mx*(1-1)), 'Linewidth', 1.0, 'DisplayName', '$$I(q)$$');
% hold on
% plot(q, IqCor{1}, 'Linestyle', '-', ... 'Color', 'black');
%     'Color', Color(3+7.5/mx*(2-1)), 'Linewidth', 1.0, 'DisplayName', '$$I_\mathrm{c}(q)$$');
%
% xlabel('$$q$$ [\AA$$^{-1}$$]');
% ylabel('Intensity [a.u.]');
% legend('show', 'Interpreter', 'latex', 'Location', 'Southeast');
% xlims = get(handles.axes6, 'xlim');
% plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility','off');
%

% figure('Name', 'Shaded Correction Degree', ...
%     'units', 'normalized', ...
%     'outerposition', [0.1 0 0.7 1], ...
%     'NumberTitle', 'off')
%
% % Plot phiq (uncorrected)
% plot(q, phiq{4}, 'Linestyle', '--', ... 'Color', 'black');
%     'Color', Color(3+7.5/mx*(1-1)), 'Linewidth', 1.0);
% hold on
%
% % Plot inverse FT of G(r<maxCorrRange) to see what is subtracted
% plot(q, f_cfM{4}, 'k', 'Linestyle', '-.', 'Linewidth', 1.0);
%
% % Plot corrected phiq
% plot(q, phiqCorM{4}, '-', ...
%     'Color', Color(3+7.5/mx*(1-1)), 'Linewidth', 1.0);
%
% xlabel('$$q$$ [\AA$$^{-1}$$]');
% ylabel('$$\phi(q)$$ [a.u.]');
% legend('show',  '$$\phi_\mathrm{u}(q)$$ (uncorrected)', ...
%                 '$$f_\mathrm{cf}(q,r_\mathrm{co})$$ (correction function)', ...
%                 '$$\phi_\mathrm{c}(q,r_\mathrm{co})=\phi_\mathrm{u}(q)-f_\mathrm{cf}(q,r_\mathrm{co})$$', ... %'\phi(q) \cdot exp(-bq^2)', ...
%         'Interpreter', 'latex', 'Location', 'Southeast', 'Fontsize', 14);
% xlims = get(handles.axes6, 'xlim');
% hold on
% plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility','off');
%
% % area(q, f_cfM{4}, 'FaceAlpha', 0.4, 'HandleVisibility', 'off', ...
% %     'FaceColor', 'red');
% area(q, phiq{4}, 'FaceAlpha', 0.4, 'HandleVisibility', 'off', ...
%     'FaceColor', 'red');
%
% % Set Fontsize
% set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

%% Plot Gr

Gr  = cell(CorrSteps,fileNumber);
GrM = cell(1,fileNumber);
GrU = cell(CorrSteps,fileNumber);  % uncorrected Gr
phiq_dampU = cell(1,fileNumber);  % uncorrected phiq_damp
for fNo = 1:fileNumber
    for iC = 1:CorrSteps
        % reduced PDF
        Gr{iC,fNo} = 8 * pi * phiq_damp{iC,fNo} * sin(q*r) * ds;
    end
    GrM{fNo} = mean(cell2mat(Gr(MeanRangeInd, fNo)));


    phiq_dampU{fNo} = phiq{fNo}'.*damping_function;
    GrU{fNo}        = 8 * pi * phiq_dampU{fNo} * sin(q*r) * ds;

end
assignin('base','Gr', Gr);
assignin('base','GrM', GrM);

cla(handles.axes7)
axes(handles.axes7);
FitAx7 = cell(1,fileNumber);
for fNo = 1:fileNumber
    FitAx7{fNo} = plot(r,GrM{fNo}, ...
        'Color', Color(3+7.5/mx*(fNo-1)));
    hold on
    FitAx7{fNo} = plot(r,GrU{fNo}, '--', ...
        'Color', Color(3+7.5/mx*(fNo-1)));
end
xlabel('$$r$$ [\AA]');
ylabel('$$G(r)$$ [\AA$$^{-2}$$]');
assignin('base','FitAx7',FitAx7);
xlim([0 rmax]);

% plot reference line at Gr=zero
xlims = get(handles.axes7, 'xlim');
plot([0 xlims(2)], [0 0],'k');
hold off

% to export results from Autofit without manual optimisation
handles.fitAll = handles.AutofitAll;
assignin('base','ResNAll', NAll);
assignin('base','Resq_fitAll', handles.q_fixAll);
guidata(hObject,handles)

%% Plot Gr and phiq evolution
%
% fNo = 4;
% figure();
%
% % Gr uncorrected
% subplot(3,1,1)
% box on;
% hold on;
% plot(r,Gr{1,fNo}, ...
%     'Color', Color(3+(1-1)/(CorrSteps-1)*9), ...
%     'Linewidth', 1.3);
%
% xlabel('$$r$$ [\AA]');
% ylabel('$$G_{\mathrm{u}}(r)$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
% xticks([0,2,4,6,8,10]);
% xlim([0 rmax]);
%
% % plot reference line at Gr=zero
% xlims = get(handles.axes7, 'xlim');
% plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');
% hold off
%
%
%
% % phiq
% subplot(3,1,2)
% for iC = 1:CorrSteps
%     % Plot phiq (uncorrected)
%     plot(q, phiqCor{iC,fNo}, 'Linestyle', '-', ...'Color', 'black');
%         'Color', Color(3+(iC-1)/(CorrSteps-1)*9), ...
%         'Linewidth', 1.3);
%     hold on
% end
%
% xlabel('$$q$$ [\AA$$^{-1}$$]');
% ylabel('$$\phi_{\mathrm{c}}(q,r_{\mathrm{co}})$$ [a.u.]', 'Interpreter', 'Latex');
% xticks([0,4,8,12,16,20]);
%
%
% % plot reference line at phiq=zero
% xlims = get(handles.axes6, 'xlim');
% hold on
% plot([0 xlims(2)], [0 0],'k', 'HandleVisibility','off');
% hold off
%
%
% % Gr
% subplot(3,1,3)
% for iC = 1:CorrSteps
%     plot(r,Gr{iC,fNo}, ...
%         'Color', Color(3+(iC-1)/(CorrSteps-1)*9), ...
%         'DisplayName', sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA', rc(iC)), ...
%         'Linewidth', 1.3);
%     hold on
% end
%
% xlabel('$$r$$ [\AA]');
% ylabel('$$G_{\mathrm{c}}(r,r_{\mathrm{co}})$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
% xticks([0,2,4,6,8,10]);
% legend('Location', 'best outside', 'Interpreter', 'Latex', 'Fontsize', 15);
% xlim([0 rmax]);
%
% % plot reference line at Gr=zero
% xlims = get(handles.axes7, 'xlim');
% plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');
% hold off



%% Export data file

filename = 'results_temp.xlsx';
delete('results_temp.xlsx');
sheetnames = {  'I(q) (raw)', 'I(q) (corrected)', ...
    'F(q) (corr.) for RMC', ...
    'F(q) norm. (corr.) for RMC', ...
    'S(q) (raw)', 'S(q) (corrected)', ...
    'G(r) (raw)', 'G(r) (corrected)', ...
    'f(q)'};

fqAll = {f1, f2, f3, f4, f5};

pre = 'f(q) for Z=';

fq_head = { [pre, num2str(elem1-1), ' (p=', num2str(e_r1), ')']; ...
    [pre, num2str(elem2-1), ' (p=', num2str(e_r2), ')']; ...
    [pre, num2str(elem3-1), ' (p=', num2str(e_r3), ')']; ...
    [pre, num2str(elem4-1), ' (p=', num2str(e_r4), ')']; ...
    [pre, num2str(elem5-1), ' (p=', num2str(e_r5), ')']};

data = {I, IqCor, FqCor, FqNCor, Sq(1,:), SqM, Gr(1,:), GrM, fqAll};

% transpose content of Sq(1,:), SqM, Gr(1,:), GrM, FqCor
for sh = 5:length(sheetnames)-1
    data{sh} = cellfun(@transpose, data{sh}, 'UniformOutput', false);
end

% create several sheets in excel file
for sh = 1:length(sheetnames)
    % check if data is in real or reciprocal space
    if isempty(strfind(sheetnames{sh},'q'))
        xdata = r';
        xhead = 'r(Ang)';
    else
        xdata = q;
        xhead = 'q(1/Ang)';
    end
    ydata = cell2mat(data{sh});

    % define first line
    if sheetnames(sh) ~= "f(q)"
        head = [xhead; TagName']';
    else
        head = [xhead, fq_head'];
    end

    % write into excel file
    writecell(head, filename, 'Sheet', sheetnames{sh}, 'Range', 'A1');
    writematrix(xdata, filename, 'Sheet', sheetnames{sh}, 'Range', 'A2');
    writematrix(ydata, filename, 'Sheet', sheetnames{sh}, 'Range', 'B2');
end
disp('results_temp.xlsx generated')


%% -----------------------------------------------------------------------
% ------------------------ CORRECT ALL ------------------------------------
% -------------------------------------------------------------------------

% ---------------------------- PLUS ---------------------------------------
function Button_Correct_All_Plus_Callback(hObject, eventdata, handles)
realSelectedFitCurve = evalin('base', 'SelectedFitCurve');
fileNumber = handles.fileNumber;
DeltaN = str2double(get(handles.edit_DeltaN, 'String'));

ResSSAll = evalin('base', 'ResSSAll');
ResNAll = evalin('base', 'ResNAll');
Resq_fitAll = evalin('base', 'Resq_fitAll');
for fNo = 1:fileNumber
    handles.NAll{fNo} = ResNAll{fNo} + DeltaN;
    handles.q_fixAll{fNo} = Resq_fitAll{fNo};
end
for fNo = 1:fileNumber
    set(handles.popSelectCurve, 'Value', fNo);
    guidata(hObject,handles);
    assignin('base','SelectedFitCurve', fNo);
    Button_Correct_Selected_Callback(handles.Button_Correct_Selected,...
        eventdata, handles);
    guidata(hObject,handles);
    ResNAll = evalin('base', 'ResNAll');
    Resq_fitAll = evalin('base', 'Resq_fitAll');
    guidata(hObject,handles);
end
% reset to previous selected curve
SelectedFitCurve = realSelectedFitCurve;
set(handles.popSelectCurve, 'Value', SelectedFitCurve);
set(handles.text_SS, 'String', ResSSAll{SelectedFitCurve});
set(handles.text_N, 'String', ResNAll{SelectedFitCurve});
set(handles.edit_N, 'String', handles.NAll{SelectedFitCurve});
set(handles.text_q_fit, 'String', Resq_fitAll{SelectedFitCurve});
set(handles.edit_q_fit, 'String', handles.q_fixAll{SelectedFitCurve});
assignin('base','SelectedFitCurve', SelectedFitCurve);

guidata(hObject,handles)

% ---------------------------- MINUS --------------------------------------
function Button_Correct_All_Minus_Callback(hObject, eventdata, handles)
realSelectedFitCurve = evalin('base', 'SelectedFitCurve');
fileNumber = handles.fileNumber;
DeltaN = str2double(get(handles.edit_DeltaN, 'String'));

ResNAll = evalin('base', 'ResNAll');
Resq_fitAll = evalin('base', 'Resq_fitAll');
for fNo = 1:fileNumber
    handles.NAll{fNo} = ResNAll{fNo} - DeltaN;
    handles.q_fixAll{fNo} = Resq_fitAll{fNo};
end
for fNo = 1:fileNumber
    set(handles.popSelectCurve, 'Value', fNo);
    guidata(hObject,handles);
    assignin('base','SelectedFitCurve', fNo);
    Button_Correct_Selected_Callback(handles.Button_Correct_Selected,...
        eventdata, handles);
    guidata(hObject,handles);
end
% reset to previous selected curve
SelectedFitCurve = realSelectedFitCurve;
ResSSAll = evalin('base', 'ResSSAll');
ResNAll = evalin('base', 'ResNAll');
Resq_fitAll = evalin('base', 'Resq_fitAll');
set(handles.popSelectCurve, 'Value', SelectedFitCurve);
set(handles.text_SS, 'String', ResSSAll{SelectedFitCurve});
set(handles.text_N, 'String', ResNAll{SelectedFitCurve});
set(handles.edit_N, 'String', handles.NAll{SelectedFitCurve});
set(handles.text_q_fit, 'String', Resq_fitAll{SelectedFitCurve});
set(handles.edit_q_fit, 'String', handles.q_fixAll{SelectedFitCurve});
assignin('base','SelectedFitCurve', SelectedFitCurve);

guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ------------------------ CORRECT SELECTED -------------------------------
% -------------------------------------------------------------------------
function Button_Correct_Selected_Callback(hObject, eventdata, handles)
%% find maximum q value smaller than edit_q_fit
SFC = evalin('base', 'SelectedFitCurve');
tri = delaunayn(handles.q);
q_index = dsearchn(handles.q,tri,handles.q_fixAll{SFC});
q_fit = handles.q(q_index);
% display q value at which fitting is done
set(handles.text_q_fit, 'String', q_fit);
set(handles.edit_q_fit, 'String', q_fit);

handles.q_fit = q_fit;

guidata(hObject,handles)

% ---------------------------------------------------------------------

%% Compute fq_sq = <f(s)>^2
s2 = handles.s2;
elem1 = handles.elem1;
elem2 = handles.elem2;
elem3 = handles.elem3;
elem4 = handles.elem4;
elem5 = handles.elem5;
e1 = handles.e1;
e2 = handles.e2;
e3 = handles.e3;
e4 = handles.e4;
e5 = handles.e5;
% compute atomic ratio (composition)
e_tot = e1+e2+e3+e4+e5;
e_r1 = e1/e_tot;
e_r2 = e2/e_tot;
e_r3 = e3/e_tot;
e_r4 = e4/e_tot;
e_r5 = e5/e_tot;
handles.e_tot = e_tot;
handles.e_r1 = e_r1;
handles.e_r2 = e_r2;
handles.e_r3 = e_r3;
handles.e_r4 = e_r4;
handles.e_r5 = e_r5;
guidata(hObject,handles)

if handles.param_val == 3 % Lobato

    paramL = handles.paramL;
    paramL_1 = paramL(elem1,:);
    paramL_2 = paramL(elem2,:);
    paramL_3 = paramL(elem3,:);
    paramL_4 = paramL(elem4,:);
    paramL_5 = paramL(elem5,:);
    paramL_elem = [paramL_1;paramL_2;paramL_3;paramL_4;paramL_5];
    paramL_table = array2table(paramL_elem,'VariableNames',{'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5'},'RowNames',{'elem1','elem2','elem3','elem4','elem5'});
    handles.paramL_table = paramL_table;

    A1_1 = paramL_table{'elem1','A1'};
    A2_1 = paramL_table{'elem1','A2'};
    A3_1 = paramL_table{'elem1','A3'};
    A4_1 = paramL_table{'elem1','A4'};
    A5_1 = paramL_table{'elem1','A5'};
    B1_1 = paramL_table{'elem1','B1'};
    B2_1 = paramL_table{'elem1','B2'};
    B3_1 = paramL_table{'elem1','B3'};
    B4_1 = paramL_table{'elem1','B4'};
    B5_1 = paramL_table{'elem1','B5'};

    f1 = ((s2*B1_1+1).^2).\(A1_1*(s2*B1_1+2))+((s2*B2_1+1).^2).\(A2_1*(s2*B2_1+2))+((s2*B3_1+1).^2).\(A3_1*(s2*B3_1+2))+((s2*B4_1+1).^2).\(A4_1*(s2*B4_1+2))+((s2*B5_1+1).^2).\(A5_1*(s2*B5_1+2));

    A1_2 = paramL_table{'elem2','A1'};
    A2_2 = paramL_table{'elem2','A2'};
    A3_2 = paramL_table{'elem2','A3'};
    A4_2 = paramL_table{'elem2','A4'};
    A5_2 = paramL_table{'elem2','A5'};
    B1_2 = paramL_table{'elem2','B1'};
    B2_2 = paramL_table{'elem2','B2'};
    B3_2 = paramL_table{'elem2','B3'};
    B4_2 = paramL_table{'elem2','B4'};
    B5_2 = paramL_table{'elem2','B5'};

    f2 = ((s2*B1_2+1).^2).\(A1_2*(s2*B1_2+2))+((s2*B2_2+1).^2).\(A2_2*(s2*B2_2+2))+((s2*B3_2+1).^2).\(A3_2*(s2*B3_2+2))+((s2*B4_2+1).^2).\(A4_2*(s2*B4_2+2))+((s2*B5_2+1).^2).\(A5_2*(s2*B5_2+2));

    A1_3 = paramL_table{'elem3','A1'};
    A2_3 = paramL_table{'elem3','A2'};
    A3_3 = paramL_table{'elem3','A3'};
    A4_3 = paramL_table{'elem3','A4'};
    A5_3 = paramL_table{'elem3','A5'};
    B1_3 = paramL_table{'elem3','B1'};
    B2_3 = paramL_table{'elem3','B2'};
    B3_3 = paramL_table{'elem3','B3'};
    B4_3 = paramL_table{'elem3','B4'};
    B5_3 = paramL_table{'elem3','B5'};

    f3 = ((s2*B1_3+1).^2).\(A1_3*(s2*B1_3+2))+((s2*B2_3+1).^2).\(A2_3*(s2*B2_3+2))+((s2*B3_3+1).^2).\(A3_3*(s2*B3_3+2))+((s2*B4_3+1).^2).\(A4_3*(s2*B4_3+2))+((s2*B5_3+1).^2).\(A5_3*(s2*B5_3+2));

    A1_4 = paramL_table{'elem4','A1'};
    A2_4 = paramL_table{'elem4','A2'};
    A3_4 = paramL_table{'elem4','A3'};
    A4_4 = paramL_table{'elem4','A4'};
    A5_4 = paramL_table{'elem4','A5'};
    B1_4 = paramL_table{'elem4','B1'};
    B2_4 = paramL_table{'elem4','B2'};
    B3_4 = paramL_table{'elem4','B3'};
    B4_4 = paramL_table{'elem4','B4'};
    B5_4 = paramL_table{'elem4','B5'};

    f4 = ((s2*B1_4+1).^2).\(A1_4*(s2*B1_4+2))+((s2*B2_4+1).^2).\(A2_4*(s2*B2_4+2))+((s2*B3_4+1).^2).\(A3_4*(s2*B3_4+2))+((s2*B4_4+1).^2).\(A4_4*(s2*B4_4+2))+((s2*B5_4+1).^2).\(A5_4*(s2*B5_4+2));

    A1_5 = paramL_table{'elem5','A1'};
    A2_5 = paramL_table{'elem5','A2'};
    A3_5 = paramL_table{'elem5','A3'};
    A4_5 = paramL_table{'elem5','A4'};
    A5_5 = paramL_table{'elem5','A5'};
    B1_5 = paramL_table{'elem5','B1'};
    B2_5 = paramL_table{'elem5','B2'};
    B3_5 = paramL_table{'elem5','B3'};
    B4_5 = paramL_table{'elem5','B4'};
    B5_5 = paramL_table{'elem5','B5'};

    f5 = ((s2*B1_5+1).^2).\(A1_5*(s2*B1_5+2))+((s2*B2_5+1).^2).\(A2_5*(s2*B2_5+2))+((s2*B3_5+1).^2).\(A3_5*(s2*B3_5+2))+((s2*B4_5+1).^2).\(A4_5*(s2*B4_5+2))+((s2*B5_5+1).^2).\(A5_5*(s2*B5_5+2));

else % Kirkland

    paramK = handles.paramK;
    paramK_1 = paramK(elem1,:);
    paramK_2 = paramK(elem2,:);
    paramK_3 = paramK(elem3,:);
    paramK_4 = paramK(elem4,:);
    paramK_5 = paramK(elem5,:);
    paramK_elem = [paramK_1;paramK_2;paramK_3;paramK_4;paramK_5];
    paramK_table = array2table(paramK_elem,'VariableNames',{'a1','b1','a2','b2','a3','b3','c1','d1','c2','d2','c3','d3'},'RowNames',{'elem1','elem2','elem3','elem4','elem5'});
    handles.paramK_table = paramK_table;

    a1_1 = paramK_table{'elem1','a1'};
    a2_1 = paramK_table{'elem1','a2'};
    a3_1 = paramK_table{'elem1','a3'};
    b1_1 = paramK_table{'elem1','b1'};
    b2_1 = paramK_table{'elem1','b2'};
    b3_1 = paramK_table{'elem1','b3'};
    c1_1 = paramK_table{'elem1','c1'};
    c2_1 = paramK_table{'elem1','c2'};
    c3_1 = paramK_table{'elem1','c3'};
    d1_1 = paramK_table{'elem1','d1'};
    d2_1 = paramK_table{'elem1','d2'};
    d3_1 = paramK_table{'elem1','d3'};

    f1 = ((s2+b1_1).\a1_1)+((s2+b2_1).\a2_1)+((s2+b3_1).\a3_1)+(exp(-s2.*d1_1).*c1_1)+(exp(-s2.*d2_1).*c2_1)+(exp(-s2.*d3_1).*c3_1);

    a1_2 = paramK_table{'elem2','a1'};
    a2_2 = paramK_table{'elem2','a2'};
    a3_2 = paramK_table{'elem2','a3'};
    b1_2 = paramK_table{'elem2','b1'};
    b2_2 = paramK_table{'elem2','b2'};
    b3_2 = paramK_table{'elem2','b3'};
    c1_2 = paramK_table{'elem2','c1'};
    c2_2 = paramK_table{'elem2','c2'};
    c3_2 = paramK_table{'elem2','c3'};
    d1_2 = paramK_table{'elem2','d1'};
    d2_2 = paramK_table{'elem2','d2'};
    d3_2 = paramK_table{'elem2','d3'};

    f2 = ((s2+b1_2).\a1_2)+((s2+b2_2).\a2_2)+((s2+b3_2).\a3_2)+(exp(-s2.*d1_2).*c1_2)+(exp(-s2.*d2_2).*c2_2)+(exp(-s2.*d3_2).*c3_2);

    a1_3 = paramK_table{'elem3','a1'};
    a2_3 = paramK_table{'elem3','a2'};
    a3_3 = paramK_table{'elem3','a3'};
    b1_3 = paramK_table{'elem3','b1'};
    b2_3 = paramK_table{'elem3','b2'};
    b3_3 = paramK_table{'elem3','b3'};
    c1_3 = paramK_table{'elem3','c1'};
    c2_3 = paramK_table{'elem3','c2'};
    c3_3 = paramK_table{'elem3','c3'};
    d1_3 = paramK_table{'elem3','d1'};
    d2_3 = paramK_table{'elem3','d2'};
    d3_3 = paramK_table{'elem3','d3'};

    f3 = ((s2+b1_3).\a1_3)+((s2+b2_3).\a2_3)+((s2+b3_3).\a3_3)+(exp(-s2.*d1_3).*c1_3)+(exp(-s2.*d2_3).*c2_3)+(exp(-s2.*d3_3).*c3_3);

    a1_4 = paramK_table{'elem4','a1'};
    a2_4 = paramK_table{'elem4','a2'};
    a3_4 = paramK_table{'elem4','a3'};
    b1_4 = paramK_table{'elem4','b1'};
    b2_4 = paramK_table{'elem4','b2'};
    b3_4 = paramK_table{'elem4','b3'};
    c1_4 = paramK_table{'elem4','c1'};
    c2_4 = paramK_table{'elem4','c2'};
    c3_4 = paramK_table{'elem4','c3'};
    d1_4 = paramK_table{'elem4','d1'};
    d2_4 = paramK_table{'elem4','d2'};
    d3_4 = paramK_table{'elem4','d3'};

    f4 = ((s2+b1_4).\a1_4)+((s2+b2_4).\a2_4)+((s2+b3_4).\a3_4)+(exp(-s2.*d1_4).*c1_4)+(exp(-s2.*d2_4).*c2_4)+(exp(-s2.*d3_4).*c3_4);

    a1_5 = paramK_table{'elem5','a1'};
    a2_5 = paramK_table{'elem5','a2'};
    a3_5 = paramK_table{'elem5','a3'};
    b1_5 = paramK_table{'elem5','b1'};
    b2_5 = paramK_table{'elem5','b2'};
    b3_5 = paramK_table{'elem5','b3'};
    c1_5 = paramK_table{'elem5','c1'};
    c2_5 = paramK_table{'elem5','c2'};
    c3_5 = paramK_table{'elem5','c3'};
    d1_5 = paramK_table{'elem5','d1'};
    d2_5 = paramK_table{'elem5','d2'};
    d3_5 = paramK_table{'elem5','d3'};

    f5 = ((s2+b1_5).\a1_5)+((s2+b2_5).\a2_5)+((s2+b3_5).\a3_5)+(exp(-s2.*d1_5).*c1_5)+(exp(-s2.*d2_5).*c2_5)+(exp(-s2.*d3_5).*c3_5);
end

fq = (f1.*e_r1) + (f2.*e_r2) + (f3.*e_r3) + (f4.*e_r4) + (f5.*e_r5);
fq_sq = fq.^2;

handles.fq_sq = fq_sq;
guidata(hObject,handles)
% -----------------------------------------------------------------------

%% Compute gq = <f^2(q)>
gq = (f1.^2*e_r1) + (f2.^2*e_r2) + (f3.^2*e_r3) + (f4.^2*e_r4) + (f5.^2*e_r5);
handles.gq = gq;
guidata(hObject,handles)
% -----------------------------------------------------------------------

% Compute fitting parameter C
f = find(handles.q == handles.q_fit);
C = handles.I{SFC}(f) - handles.gq(f)*handles.NAll{SFC};
% Compute fitting curve edit_N*gq+C
fit = handles.NAll{SFC}*handles.gq + C;
handles.fitAll{SFC} = fit;

% Update corrected fit curve plot
axes(handles.axes4);
FitAx4 = evalin('base','FitAx4');
set(FitAx4{SFC}, 'XData', handles.q, 'YData', handles.fitAll{SFC});
assignin('base','FitAx4',FitAx4);

% Update magnified view
L = uint16(length(handles.q));
q2 = handles.q(L/2:L);
handles.q2 = q2;
handles.mag_fit = handles.fitAll{SFC}(L/2:L);

axes(handles.axes5);
FitAx5 = evalin('base','FitAx5');
set(FitAx5{SFC}, 'XData', handles.q2, 'YData', handles.mag_fit);
assignin('base','FitAx5',FitAx5);

% Update goodness of fit (sum of squared deviations)
ResSSAll      = evalin('base','ResSSAll');
ResDamp       = evalin('base','ResDamp');
ResSSAll{SFC} = sprintf('%12.6e',...
    sum((handles.I{SFC}-handles.fitAll{SFC}).^2));
ResDamp{SFC}  = handles.edit_damping;
assignin('base','ResSSAll',ResSSAll);
assignin('base','ResDamp',ResDamp);
assignin('base','ResNAll', handles.NAll);
assignin('base','Resq_fitAll', handles.q_fixAll);

set(handles.text_SS, 'String', ResSSAll{SFC});
set(handles.text_N, 'String', handles.NAll{SFC});
set(handles.edit_N, 'String', handles.NAll{SFC});
set(handles.text_damping, 'String', handles.edit_damping);

guidata(hObject,handles)
% ----------------------------------------------------------------------
q       = handles.q;
I       = handles.I;
ds      = handles.ds;
fq_sq   = handles.fq_sq;
NAll    = handles.NAll;
fitAll  = handles.fitAll; % fit = N*gq+C
s       = handles.s;
s2      = handles.s2;
damping = handles.edit_damping;
% ----------------------------------------------------------------------
% Compute phiq etc.
phiq         = evalin('base', 'phiq');
CorrSteps    = evalin('base', 'CorrSteps');
numMean      = evalin('base', 'numMean');
maxCorrRange = evalin('base', 'maxCorrRange');
CStep        = maxCorrRange/CorrSteps;  % in \AA, not in r-axis units (dr)

FTphi        = evalin('base', 'FTphi');

f_cf         = evalin('base', 'f_cf');
phiqCor      = evalin('base', 'phiqCor');
Sq           = evalin('base', 'Sq');
phiq_damp    = evalin('base', 'phiq_damp');

f_cfM         = evalin('base', 'f_cfM');
phiqCorM      = evalin('base', 'phiqCorM');
SqM           = evalin('base', 'SqM');
phiq_dampM    = evalin('base', 'phiq_dampM');

IqCor         = evalin('base', 'IqCor');
FqCor         = evalin('base', 'FqCor');
FqNCor        = evalin('base', 'FqNCor');

CDegree      = evalin('base', 'CDegree');
CDegreeM     = evalin('base', 'CDegreeM');
CDegreeME    = evalin('base', 'CDegreeME');
CDegreeMErel = CDegreeME./CDegreeM;

% Update phiq (uncorrected and selected)
phiq{SFC} = ((I{SFC} - fitAll{SFC}).*s) ./ (NAll{SFC}*fq_sq);

% Correct high-frequency noise of phiq for q>qnoise (mean interpolation)
% rmaxhf = 20;
% drhf   = 0.01;
% rhf    = 0.01:drhf:rmaxhf;

qnoise = 3.;
wexp   = 2;
[~,iqnoise] = min(abs(q-qnoise));

% FTphiu   = 8*pi      *phiq{SFC}' *sin(q*rhf)  *ds;  % FT of uncorrected phiq
% iFTphiu  = 1./(2*pi) *FTphiu  *sin(q*rhf)' *drhf;  % iFT of uncorrected phiq
iFTphiu  = smooth(phiq{SFC}, 35);
phiqHfC = phiq{SFC};
for iq = iqnoise:length(q)
    phiqHfC(iq) = (phiq{SFC}(iq)+iFTphiu(iq)*(q(iq)-qnoise)^wexp)/ ...
        (1+(q(iq)-qnoise)^wexp);
end
phiq{SFC} = phiqHfC;

assignin('base','phiq', phiq);

%% Correct phiq by subtracting back-FT of G(r<maxCorrRange)
r            = evalin('base', 'r');
dr           = evalin('base', 'dr');
MeanRangeInd = evalin('base', 'MeanRangeInd');
% Forth- and back-FT of G(r<1A)
FTphi{SFC} = 8 * pi * phiq{SFC}'* sin(q*r) * ds ;
FTSys      = 0*FTphi{SFC};

% qmax: point where fitted curve crosses Iq curve
[qmax,qpos] = max(q); %value and array position of qmax

Lorch = sin(pi*q./qmax)./(pi*q./qmax);

%damping_function = exp(-s2'.*damping);
damping_function = Lorch';

% iC=1 -> no Correction, iC=CorrSteps -> max Correction
for iC = 1:CorrSteps
    CorrRange        = 1:round((((iC-1)*CStep)/dr));  % in r-axis units
    FTSys(CorrRange) = FTphi{SFC}(CorrRange);
    f_cf{iC,SFC} = 1./(2*pi) * FTSys * sin(q*r)' * dr;
    phiqCor{iC,SFC}  = phiq{SFC}' - f_cf{iC,SFC};
    Sq{iC,SFC}       = phiqCor{iC,SFC}./q'+1;

    phiq_damp{iC,SFC} = phiqCor{iC,SFC}.*damping_function;
end

% Compute averages over numMean largest correction ranges
f_cfM{SFC}  = mean(cell2mat(f_cf(MeanRangeInd, SFC)));
phiqCorM{SFC}   = mean(cell2mat(phiqCor(MeanRangeInd, SFC)));
SqM{SFC}        = mean(cell2mat(Sq(MeanRangeInd, SFC)));
phiq_dampM{SFC} = mean(cell2mat(phiq_damp(MeanRangeInd, SFC)));

% compute corrected Bragg intensity (I(q) corrected)
IqCor{SFC} = fitAll{SFC}' + (phiqCorM{SFC}*NAll{SFC}.*fq_sq')./s';
IqCor{SFC} = IqCor{SFC}';

% compute corrected F(q) and normalized F(q) for RMC fit (notation as in RMCProfile manual)
%FqCor{SFC}  = fq_sq.*(SqM{SFC}'-1);
FqCor{SFC}  = gq.*(SqM{SFC}'-1);
FqNCor{SFC} = SqM{SFC}'-1;

%% Compute Correction Degree (CDegree)
A = evalin('base', 'A');
B = evalin('base', 'B');

for iC = 1:CorrSteps
    A(iC,SFC) = trapz(q, abs(f_cf{iC,SFC}));
end
B(1,SFC) = trapz(q, abs(phiq{SFC}));

CDegree   = A./B;
CDegreeM  = mean(CDegree(MeanRangeInd,:),1);
CDegreeME = std(CDegree(MeanRangeInd,:),1);
CDegreeMErel = CDegreeME./CDegreeM;

%% Update workspace
assignin('base', 'phiq', phiq);
assignin('base', 'FTphi', FTphi);

assignin('base', 'f_cf', f_cf);
assignin('base', 'phiqCor', phiqCor);
assignin('base', 'Sq', Sq);
assignin('base', 'phiq_damp', phiq_damp);

assignin('base', 'f_cfM', f_cfM);
assignin('base', 'phiqCorM', phiqCorM);
assignin('base', 'SqM', SqM);
assignin('base', 'phiq_dampM', phiq_dampM);

assignin('base', 'IqCor', IqCor);
assignin('base', 'FqCor', FqCor);
assignin('base', 'FqNCor', FqNCor);

assignin('base', 'CDegree', CDegree);
assignin('base', 'CDegreeM', CDegreeM);
assignin('base', 'CDegreeME', CDegreeME);

% Update GUI display for CDegree (shows result of selected curve)
set(handles.text_CDegree, 'String', CDegreeM(SFC));
guidata(hObject,handles)

%% Update phiq Plot
axes(handles.axes6);
FitAx6a = evalin('base', 'FitAx6a');
FitAx6b = evalin('base', 'FitAx6b');
FitAx6c = evalin('base', 'FitAx6c');
FitAx6d = evalin('base', 'FitAx6d');
set(FitAx6a{SFC}, 'XData', q, 'YData', phiq{SFC});
set(FitAx6b{SFC}, 'XData', q, 'YData', phiqCorM{SFC});
set(FitAx6c{SFC}, 'XData', q, 'YData', phiq_dampM{SFC});
set(FitAx6d{SFC}, 'XData', q, 'YData', f_cfM{SFC});
assignin('base', 'FitAx6a', FitAx6a);
assignin('base', 'FitAx6b', FitAx6b);
assignin('base', 'FitAx6c', FitAx6c);
assignin('base', 'FitAx6d', FitAx6d);

% plot reference line at phiq=zero
xlims = get(handles.axes6, 'xlim');
hold on
plot([xlims(1) xlims(2)], [0 0],'k', 'HandleVisibility','off');
hold off

% ----------------------------------------------------------------------

%% Update Gr Plot

% Update reduced PDF
Gr  = evalin('base','Gr');
GrM = evalin('base','GrM');
for iC = 1:CorrSteps
    % reduced PDF
    Gr{iC,SFC} = 8 * pi * phiq_damp{iC,SFC} * sin(q*r) * ds;
end
GrM{SFC} = mean(cell2mat(Gr(MeanRangeInd, SFC)));
assignin('base','Gr',Gr);
assignin('base','GrM',GrM);

% Update Gr Plot
axes(handles.axes7);
FitAx7 = evalin('base', 'FitAx7');
set(FitAx7{SFC}, 'XData', r, 'YData', GrM{SFC});
assignin('base', 'FitAx7', FitAx7);

guidata(hObject,handles)

%% Export data file

TagName = evalin('base', 'TagName');

filename = 'results_temp.xlsx';
delete('results_temp.xlsx');
sheetnames = {  'I(q) (raw)', 'I(q) (corrected)', ...
    'F(q) (corr.) for RMC', ...
    'F(q) norm. (corr.) for RMC', ...
    'S(q) (raw)', 'S(q) (corrected)', ...
    'G(r) (raw)', 'G(r) (corrected)', ...
    'f(q)'};

fqAll = {f1, f2, f3, f4, f5};

pre = 'f(q) for Z=';

fq_head = { [pre, num2str(elem1-1), ' (p=', num2str(e_r1), ')']; ...
    [pre, num2str(elem2-1), ' (p=', num2str(e_r2), ')']; ...
    [pre, num2str(elem3-1), ' (p=', num2str(e_r3), ')']; ...
    [pre, num2str(elem4-1), ' (p=', num2str(e_r4), ')']; ...
    [pre, num2str(elem5-1), ' (p=', num2str(e_r5), ')']};

data = {I, IqCor, FqCor, FqNCor, Sq(1,:), SqM, Gr(1,:), GrM, fqAll};

% transpose content of Sq(1,:), SqM, Gr(1,:), GrM, FqCor
for sh = 5:length(sheetnames)-1
    data{sh} = cellfun(@transpose, data{sh}, 'UniformOutput', false);
end

% create several sheets in excel file
for sh = 1:length(sheetnames)
    % check if data is in real or reciprocal space
    if isempty(strfind(sheetnames{sh},'q'))
        xdata = r';
        xhead = 'r(Ang)';
    else
        xdata = q;
        xhead = 'q(1/Ang)';
    end
    ydata = cell2mat(data{sh});

    % define first line
    if sheetnames(sh) ~= "f(q)"
        head = [xhead; TagName']';
    else
        head = [xhead, fq_head'];
    end

    % write into excel file
    writecell(head, filename, 'Sheet', sheetnames{sh}, 'Range', 'A1');
    writematrix(xdata, filename, 'Sheet', sheetnames{sh}, 'Range', 'A2');
    writematrix(ydata, filename, 'Sheet', sheetnames{sh}, 'Range', 'B2');
end
disp('results_temp.xlsx generated')



%% ------------------------------------------------------------------------
% ------------------------ SAVE AND EXPORT --------------------------------
% -------------------------------------------------------------------------
function Export_ClickedCallback(hObject, eventdata, handles)
[file,path] = uiputfile('Results.xls','Export results to Excel as');

% retrieve raw data
q = handles.q;
I = handles.I;
fit = handles.fit;

phiq = evalin('base','phiq');
phiq_damp = evalin('base','phiq_damp');
r = evalin('base','r')';
Gr = evalin('base','Gr')';

T1 = table(q,I,fit,phiq,phiq_damp,...
    'VariableNames',{'q' 'I' 'fit' 'phiq' 'phiq_damp'});
T2 = table(r,Gr,...
    'VariableNames',{'r' 'Gr'});

% retrieve fitting parameters
ds = handles.ds;
q_fixed = handles.q_fix;
N = handles.N;
damping = handles.edit_damping;
EName1 = handles.EName1;
EName2 = handles.EName2;
EName3 = handles.EName3;
EName4 = handles.EName4;
EName5 = handles.EName5;
e_r1 = handles.e_r1;
e_r2 = handles.e_r2;
e_r3 = handles.e_r3;
e_r4 = handles.e_r4;
e_r5 = handles.e_r5;

if handles.param_val == 3
    Parameterisation = 'Lobato';
else
    Parameterisation = 'Kirkland';
end

P1 = {'Factor','ds','qmax','N','damping'};
P2 = {Parameterisation,ds,q_fixed,N,damping};
P3 = {EName1,EName2,EName3,EName4,EName5};
P4 = {e_r1,e_r2,e_r3,e_r4,e_r5};

C = [P1;P2;P3;P4];
T3 = cell2table(C);

% export to single Excel file with 2 sheets
filename = sprintf('%s\\%s',path,file);
writetable(T3,filename,'WriteVariableNames',0);
writetable(T1,filename,'Sheet',2);
writetable(T2,filename,'Sheet',2,'Range','F1');
% -------------------------------------------------------------------------
function Export_txt_ClickedCallback(hObject, eventdata, handles)
% --- export raw data of plots and fitting parameters to separate csv files
% --- (No Excel support in Mac)

% Ask user to choose directory to export csv files
folder = uigetdir('','Select folder to export CSV results to');

% retrieve raw data
q = handles.q;
I = handles.I;
fit = handles.fit;

phiq = evalin('base','phiq');
phiq_damp = evalin('base','phiq_damp');

r = evalin('base','r')';
Gr = evalin('base','Gr')';

T1 = table(q,I,fit,phiq,phiq_damp,...
    'VariableNames',{'q' 'I' 'fit' 'phiq' 'phiq_damp'});
T2 = table(r,Gr,...
    'VariableNames',{'r' 'Gr'});

% retrieve fitting parameters
ds = handles.ds;
q_fixed = handles.q_fix;
N = handles.N;
damping = handles.edit_damping;
EName1 = handles.EName1;
EName2 = handles.EName2;
EName3 = handles.EName3;
EName4 = handles.EName4;
EName5 = handles.EName5;
e_r1 = handles.e_r1;
e_r2 = handles.e_r2;
e_r3 = handles.e_r3;
e_r4 = handles.e_r4;
e_r5 = handles.e_r5;

if handles.param_val == 3
    Parameterisation = 'Lobato';
else
    Parameterisation = 'Kirkland';
end

P1 = {'Factor','ds','qmax','N','damping'};
P2 = {Parameterisation,ds,q_fixed,N,damping};
P3 = {EName1,EName2,EName3,EName4,EName5};
P4 = {e_r1,e_r2,e_r3,e_r4,e_r5};

C = [P1;P2;P3;P4];
T3 = cell2table(C);

% export to separate csv files
filename1 = sprintf('%s\\Results_q.csv',folder);
writetable(T1,filename1);
filename2 = sprintf('%s\\Results_r.csv',folder);
writetable(T2,filename2);
filename3 = sprintf('%s\\Parameters.csv',folder);
writetable(T3,filename3,'WriteVariableNames',0);
% -------------------------------------------------------------------------
function Save_ClickedCallback(hObject, eventdata, handles)
% Ask user to choose directory to save Iq, Phiq and Gr plots as png
folder = uigetdir('','Select folder to save plots');

set(groot,'defaultFigurePaperPositionMode','auto');

% Iq
ax1 = handles.axes4;
ax1.Units = 'pixels';
pos = ax1.Position;
ti = ax1.TightInset;
rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
PlotI = getframe(ax1,rect);

figure('visible','off','NumberTitle', 'off');
imshow(PlotI.cdata);
filename1 = sprintf('%s\\Plot_Iq',folder);
print(gcf,filename1,'-dpng');
close(gcf);

% Phiq
ax2 = handles.axes6;
ax2.Units = 'pixels';
pos = ax2.Position;
ti = ax2.TightInset;
rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
PlotPhi = getframe(ax2,rect);

figure('visible','off','NumberTitle', 'off');
imshow(PlotPhi.cdata);
filename2 = sprintf('%s\\Plot_Phiq',folder);
print(gcf,filename2,'-dpng');
close(gcf);

% Gr
ax3 = handles.axes7;
ax3.Units = 'pixels';
pos = ax3.Position;
ti = ax3.TightInset;
rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
PlotGr = getframe(ax3,rect);

figure('visible','off','NumberTitle', 'off');
imshow(PlotGr.cdata);
filename3 = sprintf('%s\\Plot_Gr',folder);
print(gcf,filename3,'-dpng');
close(gcf);

%% ------------------------------------------------------------------------
% -------------------- SAVE AND EXPORT OPTIONS ----------------------------
% -------------------------------------------------------------------------
function checkbox_save_Callback(hObject, eventdata, handles)
handles.save = get(hObject,'Value');
guidata(hObject,handles)
function checkbox_overwrite_Callback(hObject, eventdata, handles)
handles.overwrite = get(hObject,'Value');
guidata(hObject,handles)



%% ------------------------------------------------------------------------
% ---------------------------- ANALYSE S(q) -------------------------------
% -------------------------------------------------------------------------
function push_Sq2Sq1_Callback(hObject, eventdata, handles)
fileNumber = handles.fileNumber;
foldername = handles.datfname;
TagName    = evalin('base', 'TagName');
TagNameLegend = evalin('base', 'TagNameLegend');
TagTitle   = evalin('base', 'TagTitle');
Sq         = evalin('base', 'Sq');
SqM        = evalin('base', 'SqM');
CorrSteps  = evalin('base', 'CorrSteps');
rc         = evalin('base', 'rc');
CDegree    = evalin('base', 'CDegree');
numMean    = evalin('base', 'numMean');  % number of correction ranges over which is averaged

rco          = evalin('base', 'rco');
rcoind       = evalin('base', 'rcoind');
MeanRangeInd = evalin('base', 'MeanRangeInd');

q  = handles.q;
mx = handles.mx;


% -------------------------------------------------------------------------
%% S(q) & Delta S(q) = S(q)-S_0(q)

% Calculate Delta S(q)
ref = 1;
DSq = cell(1,fileNumber);
for fNo = 1:fileNumber
    DSq{fNo} = SqM{fNo}-SqM{ref};
end

% Initialize Plot
SqPlot = figure('Name', 'S(q)', 'units', 'normalized', ...
    'NumberTitle', 'off', 'outerposition', [0.1 0 0.7 1]);

% Plot S(q)
subplot(2,1,1);
%subplot(1,1,1);
grid on
hold on
box on
for fNo = 1:fileNumber
    if fNo == ref
        plot(q, SqM{fNo}, 'Linestyle', '-',...
            'DisplayName', [TagNameLegend{fNo}(:)', ' $$(S_0)$$'],...
            'Color', Color(3+7.5/mx*(fNo-1)), ...
            'Linewidth', 1.1);
    else
        plot(q, SqM{fNo}, 'Linestyle', '-',...
            'DisplayName', TagNameLegend{fNo}(:)',...
            'Color', Color(3+7.5/mx*(fNo-1)), ...
            'Linewidth', 1.1);
    end
end
hold off
xlabel('$q$ [\AA$^{-1}]$', 'Interpreter', 'Latex');
ylabel('$S(q)$ [a.u.]', 'Interpreter', 'Latex');
legend('show', 'Location', 'NorthEast', 'Interpreter', 'Latex');

% Plot reference line at Sq=1
xlims = get(handles.axes6, 'xlim');
hold on
plot([0 xlims(2)], [1 1], 'k', 'HandleVisibility','off');
hold off

% Plot Delta S(q)
subplot(2,1,2);
grid on
hold on
box on
for fNo = 1:fileNumber
    plot(q, DSq{fNo}, ...
        'DisplayName', TagNameLegend{fNo}(:)',...
        'Color', Color(3+7.5/mx*(fNo-1)), ...
        'Linewidth', 1.2);
end
hold off
xlabel('$q$ [\AA$^{-1}]$', 'Interpreter', 'Latex');
ylabel('$\Delta S_i(q) = S_i(q) - S_0(q)$', 'Interpreter', 'Latex');

% plot reference line at Delta Sq=0
% xlims = get(handles.axes6, 'xlim');
% ylims = get(gca, 'ylim');
% ylim(1.4*ylims);
hold on
plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility','off');
hold off

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 16)

%% Fit first two peaks
numPeaks = 2;  % peak = maximum

mu       = cell(1,fileNumber);
sigma    = cell(1,fileNumber);
scaling  = cell(1,fileNumber);
muE      = cell(1,fileNumber);
sigmaE   = cell(1,fileNumber);
scalingE = cell(1,fileNumber);
SqScaling = cell(1,fileNumber);

FitSqPlot = cell(fileNumber);
SqEvoPlot = cell(fileNumber);
% Computing fit and plotting evolution of S(q)
for fNo = 1:fileNumber
    Sq0 = cellfun(@(x) x-1, Sq(:,fNo), 'UniformOutput', 0);  % oscillates about 0 for fit

    SqEvoPlot{fNo} = figure('Name', 'S(q) Evolution upon correction ' ...
        + TagTitle + ' ' + TagName{fNo}, ...
        'units', 'normalized', ...
        'outerposition', [0.1 0 0.7 1], ...
        'NumberTitle', 'off');
    title(TagNameLegend{fNo}, 'Interpreter', 'Latex');
    xlabel('$$q$$ [\AA$$^{-1}$$]', 'Interpreter', 'Latex');
    ylabel('$$S(q,r_{\mathrm{co}})$$', 'Interpreter', 'Latex');
    hold on
    box on

    Sq01 = cellfun(@(x) x, Sq(:,fNo), 'UniformOutput', 0);
    for iC = 1:CorrSteps
        plot(q, Sq01{iC}, 'Color', Color(3+7.5/(CorrSteps-1)*(iC-1)), ...
            'DisplayName', ...
            ... sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA$$\\,\\,$$($$C_\\mathrm{D}=%.2f$$)', ...
            ... rc(iC), CDegree(iC,fNo)), ...
            sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA', ...
            rc(iC)), ...
            'Linewidth', 1.0);
    end


    % Plot reference line at Sq=0
    xlims = get(handles.axes6, 'xlim');
    hold on
    plot([0 xlims(2)], [1 1], 'k', 'HandleVisibility', 'off');
    hold off
    legend('Interpreter', 'Latex', 'Location', 'southeast')

    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

    assignin('base', 'Sq01', Sq01);
    assignin('base', 'q', q);


    [   mu{fNo},  sigma{fNo},  scaling{fNo}, ...
        ... error on individual- (and not mean-) value: i.e. not devided by sqrt(number of different fitdata lengths around maximum)
        muE{fNo}, sigmaE{fNo}, scalingE{fNo}, ...
        FitSqPlot{fNo}] = ...
        Gaussfit_Sq( q', Sq0, numPeaks, TagNameLegend{fNo}(:)' );
    scaling{fNo}   = scaling{fNo} - 0.05;  % Gaussfit baseline at 0.95
    SqScaling{fNo} = scaling{fNo} + 1;
end

% Compute corrected averages of numMean correction ranges
% mu--> 1x3-cell 15x2-double
muM       = cell(1,fileNumber);  % muM-> 1x3-cell  1x2-double
sigmaM    = cell(1,fileNumber);
scalingM  = cell(1,fileNumber);
muEM      = cell(1,fileNumber);
sigmaEM   = cell(1,fileNumber);
scalingEM = cell(1,fileNumber);
SqScalingM = cell(1,fileNumber);
for fNo = 1:fileNumber
    muM{fNo}       = mean(mu{fNo}(MeanRangeInd, :));
    muEM{fNo}      = std(mu{fNo}(MeanRangeInd, :));
    sigmaM{fNo}    = mean(sigma{fNo}(MeanRangeInd, :));
    sigmaEM{fNo}   = std(sigma{fNo}(MeanRangeInd, :));
    scalingM{fNo}  = mean(scaling{fNo}(MeanRangeInd, :));
    scalingEM{fNo} = std(scaling{fNo}(MeanRangeInd, :));
    SqScalingM{fNo} = mean(SqScaling{fNo}(MeanRangeInd, :));
end

%% Define plot colors
cmu     = [0.00,0.45,0.74];   % green color for mu
csc     = Color(12.375);  % purpe color for scaling
csi     = Color(4.875);   % blue color for sigma

%% Plot correction evolution of Mu/Sigma/Scaling
PeakCorrectionEvoPlot = cell(fileNumber,1);
for fNo = 1:fileNumber
    PeakCorrectionEvoPlot{fNo} = figure('Name', ...
        'S(q) Mu/Sigma/Scaling Evolution upon correction: ' ...
        + TagTitle + ' ' + TagName{fNo}, ...
        'units', 'normalized', ...
        'outerposition', [0.1 0 0.7 1], ...
        'NumberTitle', 'off');
    sgtitle(TagNameLegend{fNo}, 'Interpreter', 'Latex');

    for p = 1:numPeaks
        subplot(3,numPeaks,p)
        set(gca, 'Fontsize', 14);
        hold on
        box on

        title(['Maximum ', num2str(p)]);
        errorbar(mu{fNo}(:,p), muE{fNo}(:,p), 'Marker', 'x', ...
            'DisplayName', '', 'Linewidth', 0.8, ...
            'Color', cmu)
        plot_hline(muM{fNo}(p), muEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        xticklabels('')
        ylabel(sprintf('$$q_{%i}$$ [\\AA$$^{-1}$$]', p), ...
            'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:2:end))
        set(gca,'XTick',1:2:CorrSteps)
        grid on

        subplot(3,numPeaks,p+numPeaks)
        set(gca, 'Fontsize', 14);
        hold on
        box on
        errorbar(sigma{fNo}(:,p), sigmaE{fNo}(:,p), 'Marker', 'x', ...
            'DisplayName', '', 'Linewidth', 0.8, ...
            'Color', csi)
        plot_hline(sigmaM{fNo}(p), sigmaEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        xticklabels({''})
        ylabel(sprintf('$$\\sigma_{%i}$$ [\\AA$$^{-1}$$]', p), 'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:2:end))
        set(gca,'XTick',1:2:CorrSteps)
        grid on

        Ax = subplot(3,numPeaks,p+2*numPeaks);
        set(gca, 'Fontsize', 14);
        hold on
        box on
        errorbar(SqScaling{fNo}(:,p), scalingE{fNo}(:,p), 'Marker', 'x', ...
            'DisplayName', '', 'Linewidth', 0.8, ...
            'Color', csc)
        plot_hline(SqScalingM{fNo}(p), scalingEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        Ax.TickLabelInterpreter = 'none';
        ylabel(sprintf('$$S(q_{%i})$$', p), 'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:2:end))
        set(gca,'XTick',1:2:CorrSteps)
        grid on
    end
end

%% Plot correction evolution of q1/qx

q1qxCorrectionEvoPlot = figure('Name', ...
    'q1/qx Evolution upon correction: ' ...
    + TagTitle + ' ' + TagName{fNo}, ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');

% sgtitle('$q_1^\mathrm{ref}/q_1^x$ evolution upon correction', ...
%     'Interpreter', 'Latex');

%tiledlayout(fileNumber,1);
tiledlayout(1,1);
yl = [1 1];
for fNo = 1:fileNumber
    q1qx   = mu{1}(:,1) ./ mu{fNo}(:,1);  % ~V
    q1qxE  = sqrt( (muE{1}(:,1)  ./ mu{1}(:,1)).^2 ...
        + (muE{fNo}(:,1)./ mu{fNo}(:,1)).^2 ) .*abs(q1qx);
    q1qxM  = mean( q1qx(MeanRangeInd) );
    q1qxEM =  std( q1qx(MeanRangeInd) );

    if fNo == 1
        q1qxE(:) = 0;  % accounting for correlation when dividing by itself
    end

    %nexttile
    set(gca, 'FontSize', 14);
    hold on
    box on
    errorbar(q1qx, q1qxE, '-', 'Marker', 'x', ...
        'DisplayName', TagNameLegend{fNo}(:), ...
        'Color', Color(3+7.5/mx*(fNo-1)), ...
        'LineWidth', 1.0);
    plot_hline(q1qxM, q1qxEM, ...
        rcoind-numMean-0.1, rcoind+numMean+0.1, ...
        Color(3+7.5/mx*(fNo-1)))


    ylabel('$q_1^\mathrm{ref}/q_1^x$', 'Interpreter', 'Latex')
    xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
    xlim([0.5 CorrSteps+0.5])
    set(gca, 'XTickLabel', rc)
    set(gca, 'XTick', 1:CorrSteps)
    grid on
    legend('Location', 'NorthWest', 'box', 'on', 'Interpreter', 'Latex');

    % zoom to relevant y-limits
    %     if fNo > 1
    %         yVar = max(q1qxEM,max(q1qxE(MeanRangeInd, :)));
    %         ylim([q1qxM-5*yVar q1qxM+5*yVar])
    %     end

    % same y limits
    %     cylim = get(gca, 'ylim');
    %     if cylim(1) < yl(1) && cylim(1) ~= 0
    %         yl(1) = cylim(1);
    %     end
    %     if cylim(2) > yl(2) && cylim(2) ~= 2
    %         yl(2) = cylim(2);
    %     end
    % end
    % for fNo = 1:fileNumber
    %     nexttile(fNo);
    %     ylim([yl(1) yl(2)]);
end
% -------------------------------------------------------------------------

%% Plot correction evolution of Sq2/Sq1

Sq2Sq1CorrectionEvoPlot = figure('Name', ...
    'S(q2)/S(q1) Evolution upon correction: ' ...
    + TagTitle + ' ' + TagName{fNo}, ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');

% sgtitle('$$S(q_1)/S(q_2)$$ evolution upon correction', ...
%     'Interpreter', 'Latex');
% sgtitle('Relative Change in $$S(q_1)/S(q_2)$$ upon Correction', ...
%     'Interpreter', 'Latex');

%tiledlayout(fileNumber,1);
%tiledlayout(1,1);
yl = [1 1];
for fNo = 1:fileNumber
    Sq2Sq1   = SqScaling{fNo}(:,2) ./ SqScaling{fNo}(:,1);
    Sq2Sq1E  = sqrt( (scalingE{fNo}(:,1)./ SqScaling{fNo}(:,1)).^2 ...
        + (scalingE{fNo}(:,2)./ SqScaling{fNo}(:,2)).^2 ) ...
        .* abs(Sq2Sq1);  % Fitting errors
    Sq2Sq1M  = mean( Sq2Sq1(MeanRangeInd) );
    Sq2Sq1EM =  std( Sq2Sq1(MeanRangeInd) );  % Correction error

    Sq2Sq1rel   = (Sq2Sq1./Sq2Sq1(rcoind)-1)*100;
    Sq2Sq1Erel  = (Sq2Sq1E./Sq2Sq1)*100;
    Sq2Sq1Mrel  = (Sq2Sq1M./Sq2Sq1(rcoind)-1)*100;
    Sq2Sq1EMrel = (Sq2Sq1EM./Sq2Sq1M)*100;

    %nexttile
    set(gca, 'FontSize', 14);
    box on
    hold on


    errorbar(Sq2Sq1, Sq2Sq1E, '-', 'Marker', 'x', ...
        'DisplayName', TagNameLegend{fNo}(:), ...
        'Color', Color(3+7.5/mx*(fNo-1)), ...
        'LineWidth', 1.0);
    plot_hline(Sq2Sq1M, Sq2Sq1EM, ...
        rcoind-numMean-0.1, rcoind+numMean+0.1, ...
        Color(3+7.5/mx*(fNo-1)))
    ylabel('$S(q_2)/S(q_1)$', 'Interpreter', 'Latex')

    %     errorbar(Sq2Sq1rel,Sq2Sq1Erel, '-', 'Marker', 'x', ...
    %         'DisplayName', TagNameLegend{fNo}, ...
    %         'LineWidth', 1.0);
    %     plot_hline(Sq2Sq1Mrel, Sq2Sq1EMrel, ...
    %             rcoind-numMean-0.1, rcoind+numMean+0.1, ...
    %             Color(3+7.5/mx*(fNo-1)))
    %     ylabel('$\Delta (S_2/S_1)$ [$\%$]', 'Interpreter', 'Latex')

    %set(get(gca,'YLabel'), 'rotation', 0)
    xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
    xlim([0.5 CorrSteps+0.5])
    set(gca,'XTickLabel',rc)
    set(gca,'XTick',1:CorrSteps)
    grid on
    legend('Location', 'NorthWest', 'box', 'on', 'Interpreter', 'Latex');

    % same y limits
    cylim = get(gca, 'ylim');
    if cylim(1) < yl(1)
        yl(1) = cylim(1);
    end
    if cylim(2) > yl(2)
        yl(2) = cylim(2);
    end
end

% % scale to same y-axes
% for fNo = 1:fileNumber
%     nexttile(fNo);
%     ylim([yl(1) yl(2)]);
% end
% -------------------------------------------------------------------------

%% Plot Correction Degree (CDegree)
CDegree = evalin('base', 'CDegree');
% A = evalin('base', 'A');
% B = evalin('base', 'B');
CDegreePlot = figure('Name', ...
    'Correction Degree: Evolution upon Correction: ' + TagTitle, ...
    'units', 'normalized', ...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');
% subplot(3,1,1)
set(gca, 'Fontsize', 14);
%grid on
ylimmax = 0;
for fNo = 1:fileNumber
    %     subplot(3,1,1)
    hold on
    box on
    plot(abs(CDegree(:,fNo)), '-o', 'DisplayName', TagNameLegend{fNo}(:)', ...
        'Color', Color(3+7.5/mx*(fNo-1)), ...
        'Linewidth', 1.0);
    ylabel('$$C_\mathrm{D}$$ [a.u.]', 'Interpreter', 'Latex');
    %     subplot(3,1,2)
    %     hold on
    %     box on
    %     plot(abs(A(:,fNo)), '-o', 'DisplayName', TagNameLegend{fNo}(:)', ...
    %         'Color', Color(3+7.5/mx*(fNo-1)), ...
    %         'Linewidth', 1.0);
    %     ylabel('$$A$$ (a.u.)', 'Interpreter', 'Latex');
    %     subplot(3,1,3)
    %     hold on
    %     box on
    %     plot(abs(B(:,fNo)), '-o', 'DisplayName', TagNameLegend{fNo}(:)', ...
    %         'Color', Color(3+7.5/mx*(fNo-1)), ...
    %         'Linewidth', 1.0);
    %     ylabel('$$B$$ (a.u.)', 'Interpreter', 'Latex');
    %     Bmax    = max(abs(B(:,fNo)));
    %     ylimmax = max(ylimmax, Bmax);
    %     ylim([0,ylimmax*1.1]);
end



% plot cutoff line for r_co
yl = ylim;
plot(rcoind.*[1 1],[yl(1) yl(2)], ...
    'k', 'HandleVisibility', 'off');
plot((rcoind-numMean-0.1).*[1 1],[yl(1) yl(2)], ...
    '--k', 'HandleVisibility', 'off');
plot((rcoind+numMean+0.1).*[1 1],[yl(1) yl(2)], ...
    '--k', 'HandleVisibility', 'off');

set(gca, 'ylim', yl);
xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex');
xlim([0.5 CorrSteps+0.5])
set(gca,'XTickLabel',rc)
set(gca,'XTick',1:CorrSteps)
% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)


legend('Interpreter', 'Latex', 'Location', 'NorthWest', 'box', 'off');

%% Error propagation
for fNo = 1:fileNumber
    for p = 1:numPeaks
        muCorrErr    = muEM{fNo}(p);
        muFitErr     = max(muE{fNo}(MeanRangeInd,p));
        muEM{fNo}(p) = sqrt( muCorrErr^2 + muFitErr^2 );

        sigmaCorrErr    = sigmaEM{fNo}(p);
        sigmaFitErr     = max(sigmaE{fNo}(MeanRangeInd,p));
        sigmaEM{fNo}(p) = sqrt( sigmaCorrErr^2 + sigmaFitErr^2 );

        scalingCorrErr    = scalingEM{fNo}(p);
        scalingFitErr     = max(scalingE{fNo}(MeanRangeInd,p));
        scalingEM{fNo}(p) = sqrt( scalingCorrErr^2 + scalingFitErr^2 );
    end
end

%% Re-organize results for plotting
muM        = cell2mat(muM');
sigmaM     = cell2mat(sigmaM');
scalingM   = cell2mat(scalingM');
SqScalingM = scalingM + 1;
muEM       = cell2mat(muEM');
sigmaEM    = cell2mat(sigmaEM');
scalingEM  = cell2mat(scalingEM');

%% Plotting mu, (sigma), scaling
PeakPlot  = cell(fileNumber,1);
for p = 1:numPeaks

    PeakPlot{p} = figure('Name', ...
        'S(q) Mu/(Sigma/)Scaling: ' ...
        + TagTitle + [': Local maximum ', int2str(p)], ...
        'units','normalized', ...
        'outerposition',[0.1 0 0.7 1], ...
        'NumberTitle', 'off');

    % mu
    subplot(2,1,1)
    sgtitle(['$$S(q)$$ Maximum ', int2str(p)], ...
        'Interpreter', 'Latex');
    errorbar(muM(:,p), muEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'linewidth', 1.0, ...
        'Color', cmu)
    xticklabels('')
    ylabel(sprintf('$$q_{%i}$$ [\\AA$$^{-1}$$]', p), 'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(gca,'XTick',1:numel(TagNameLegend))
    grid on

    % sigma
    %     subplot(3,1,2)
    %     errorbar(sigmaM(:,p), sigmaEM(:,p), 'Marker', 'x', ...
    %        'DisplayName', '', 'linewidth', 1.0, ...
    %        'Color', csi)
    %     xticklabels({''})
    %     ylabel(sprintf('$$\\sigma_{%i}$$ [\\AA$$^{-1}$$]', p), 'Interpreter', 'Latex')
    %     xlim([0.5 fileNumber+0.5])
    %     set(gca,'XTick',1:numel(TagName))
    %     grid on

    % scaling
    Ax = subplot(2,1,2);
    errorbar(SqScalingM(:,p), scalingEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'linewidth', 1.0, ...
        'Color', csc)
    Ax.TickLabelInterpreter = 'none';
    ylabel(sprintf('$$S(q_{%i})$$ [a.u.]', p), 'Interpreter', 'Latex')

    xlim([0.5 fileNumber+0.5])
    set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'XTickLabel',TagNameLegend)
    set(gca,'XTick',1:numel(TagNameLegend))
    xtickangle(8)
    grid on

    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

end

%% Plotting q1, q1/qx

% Plot
q1qxPlot = figure('Name', 'S(q) Volume Change: ' + TagTitle,...
    'units', 'normalized', ...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');

subplot(2,1,1)
%sgtitle(TagTitle + ': $(q_1^\mathrm{ref}/q_1^x)^3\sim V$', 'Interpreter', 'Latex');
errorbar(muM(:,1), muEM(:,1), 'Marker', 'x', ...
    'DisplayName', '', 'linewidth', 1.0)
xticklabels({''})
ylabel('$$q_1$$ [\AA $$^{-1}$$]', 'Interpreter', 'Latex')
set(gca,'XTick',1:numel(TagName))
xlim([0.5 fileNumber+0.5])
grid on

% subplot(3,1,2)
% errorbar(muM(:,2), muEM(:,2), 'DisplayName', '')
% xticklabels({''})
% ylabel('$q_2$ [\AA $^{-1}$]', 'Interpreter', 'Latex')
% set(gca,'XTick',1:numel(TagName))
% xlim([0.5 fileNumber+0.5])
% grid on


% % q2/q1
% q2q1  = muM(:,2) ./ muM(:,1);
% q2q1E = sqrt( (muEM(:,2)./ muM(:,2)).^2 ...
%             + (muEM(:,1)./ muM(:,1)).^2 ).*abs(q2q1);
% subplot(4,1,3);
% errorbar(q2q1, q2q1E, 'DisplayName', '')
% xticklabels({''})
% ylabel('$$q_2/q_1$$', 'Interpreter', 'Latex')
% set(gca,'XTick',1:numel(TagName))
% xlim([0.5 fileNumber+0.5])
% grid on


% q1_asdep/q1_x density
q1qx  = muM(1,1) ./ muM(:,1);  % ~V
q1qxE = sqrt( (muEM(1,1)./ muM(1,1)).^2 ...
    + (muEM(:,1)./ muM(:,1)).^2 ).*abs(q1qx);
q1qxE(1) = 0;  % accounting for correlation when dividing by itself

% q1qx^3 ~ V
q1qx3  = q1qx.^3;
q1qx3E = 3*q1qx.^2.*q1qxE;

Ax = subplot(2,1,2);
errorbar(q1qx3, q1qx3E, 'Marker', 'x', ...
    'DisplayName', '', 'linewidth', 1.0)
Ax.TickLabelInterpreter = 'none';
ylabel('$(q_1^\mathrm{ref}/q_1^x)^3\sim V$', 'Interpreter', 'Latex')
set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
set(gca,'XTickLabel',TagNameLegend)
set(gca,'XTick',1:numel(TagNameLegend))
xtickangle(8)
xlim([0.5 fileNumber+0.5])
grid on

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

% -------------------------------------------------------------------------

%% Plotting S(q1), S(q2), S(q2)/S(q1)

Sq2Sq1Plot = figure('Name', 'S(q) Peak Heights ' + TagTitle,...
    'units', 'normalized', ...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');

% Plot Sq1
subplot(3,1,1)
%sgtitle(TagTitle + ': $S(q_2)/S(q_1)$', 'Interpreter', 'Latex', 'Fontsize', 15);
errorbar(SqScalingM(:,1), scalingEM(:,1), 'Marker', 'x', ...
    'DisplayName', '', 'linewidth', 1.0);
xticklabels({''})
ylabel('$$S(q_1)$$', 'Interpreter', 'Latex', 'FontSize', 14)
set(gca,'XTick',1:numel(TagName))
xlim([0.5 fileNumber+0.5])
grid on

% Plot Sq2
subplot(3,1,2)
errorbar(SqScalingM(:,2), scalingEM(:,2), 'Marker', 'x', ...
    'DisplayName', '', 'linewidth', 1.0);
xticklabels({''})
ylabel('$$S(q_2)$$', 'Interpreter', 'Latex', 'FontSize', 14)
set(gca,'XTick',1:numel(TagName))
xlim([0.5 fileNumber+0.5])
grid on

% Calculate Sq1/Sq2
Sq2Sq1  = SqScalingM(:,2) ./ SqScalingM(:,1);
Sq2Sq1E = sqrt( (scalingEM(:,2) ./ SqScalingM(:,2)).^2 ...
    + (scalingEM(:,1) ./ SqScalingM(:,1)).^2 ) .* abs(Sq2Sq1);

% Plot Sq1/Sq2
Ax = subplot(3,1,3);
errorbar(Sq2Sq1, Sq2Sq1E, 'Marker', 'x', ...
    'DisplayName', '', 'linewidth', 1.0);
Ax.TickLabelInterpreter = 'none';
ylabel('$$S(q_2)/S(q_1)$$', 'Interpreter', 'Latex', 'FontSize', 14)
set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
set(gca,'XTickLabel',TagNameLegend)
set(gca,'XTick',1:numel(TagNameLegend))
xtickangle(8)
xlim([0.5 fileNumber+0.5])
grid on

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

guidata(hObject,handles)
% -------------------------------------------------------------------------

%% save figures if checkbox "save" is activated
if handles.save == 1
    disp('saving...')
    extensions = {'fig', 'svg', 'png'};
    rehash
    % create folder name (e.g. _02_03_04_05_06_Analyse_S(q))
    res_folder_name = '_';
    for fNo = 1:fileNumber
        whitespace = find(isspace(TagName{fNo})==1, 1);
        res_folder_name = [res_folder_name, ...
            char(TagName{fNo}(1:whitespace-1)), ' '];
    end
    res_folder_name = [res_folder_name, 'Analyse_S(q)'];
    res_folder_name = strrep(res_folder_name, ' ', '_');

    % delete latest data if checkbox "overwrite" is activated
    % and if existing
    rehash
    if handles.overwrite == 1 && ...
            isfolder([handles.datpath, '\', res_folder_name]) == 1

        % find latest figure in folder
        FileInfo = dir([handles.datpath, '\', res_folder_name, '\*.fig']);
        TimeStamp = 0;
        delind = 0;
        for nfile = 1:size(FileInfo)
            FI = FileInfo(nfile,1);
            if FI.datenum > TimeStamp
                TimeStamp = FI.datenum;
                delind = nfile;
            end
        end

        % find all folders and files with same date
        FI = FileInfo(delind,1);
        sim_name = FI.name(1:20);
        for k = 1:length(extensions)
            delete([FI.folder, '\', sim_name, '*.', extensions{k}]);
        end
        delfolders = dir([FI.folder, '\', sim_name, '*']);
        for i = 1:size(delfolders)
            deldata = delfolders(i);
            delpath = [deldata.folder, '\', deldata.name];
            rmdir(delpath, 's');
        end


        % create dir if not existing
    elseif isfolder([handles.datpath, '\', res_folder_name]) == 0
        mkdir([handles.datpath, '\', res_folder_name])
    end

    % get timestamp
    date = strrep(datestr(now), ':', '-');
    date = strrep(date, ' ', '_');


    % save (SqPlot) figure (02-Apr-2020_16-36-26_Sq.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_Sq'];
    for k = 1:length(extensions)
        saveas(SqPlot, savename, extensions{k})
    end
    % save (SqEvoPlot) figure (02-Apr-2020_16-36-26_Sq_Correction_Evolution.fig)
    for fNo = 1:fileNumber
        savename = [handles.datpath, '\', res_folder_name, '\', ...
            date, '_Sq_Correction_Evolution_', strrep(TagName{fNo}, '.','')];
        for k = 1:length(extensions)
            saveas(SqEvoPlot{fNo}, savename, extensions{k})
        end
    end


    % save q1qxCorrectionEvoPlot (02-Apr-2020_16-36-26_q1qx_Correction_Evolution.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_q1qx_Correction_Evolution'];
    for k = 1:length(extensions)
        saveas(q1qxCorrectionEvoPlot, savename, extensions{k})
    end
    % save q1qxPlot (02-Apr-2020_16-36-26_q1qx.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_q1qx'];
    for k = 1:length(extensions)
        saveas(q1qxPlot, savename, extensions{k})
    end


    % save Sq2Sq1CorrectionEvoPlot (02-Apr-2020_16-36-26_Sq2Sq1_Correction_Evolution.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_Sq2Sq1_Correction_Evolution'];
    for k = 1:length(extensions)
        saveas(Sq2Sq1CorrectionEvoPlot, savename, extensions{k})
    end
    % save Sq2Sq1Plot (02-Apr-2020_16-36-26_Sq2Sq1.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_Sq2Sq1'];
    for k = 1:length(extensions)
        saveas(Sq2Sq1Plot, savename, extensions{k})
    end


    % save CDegreePlot (02-Apr-2020_16-36-26_Correction_Degree.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_Correction_Degree'];
    for k = 1:length(extensions)
        saveas(CDegreePlot, savename, extensions{k})
    end



    % create subfolder for all Peakfit results
    res_subfolder_name = [date, '_Peakfit_Results'];
    res_subfolder_path = [handles.datpath, '\', res_folder_name, '\', ...
        res_subfolder_name];
    if isfolder(res_subfolder_path) == 0
        mkdir(res_subfolder_path)
    end
    % save (PeakCorrectionEvoPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\Evolution_Peak_1+2.fig)
    for fNo = 1:fileNumber
        savename = [res_subfolder_path, '\', date, '_', ...
            strrep( strrep(TagName{fNo},' ','_'), '.',''), ...
            '_Evolution_Peak_1+2'];
        for k = 1:length(extensions)
            saveas(PeakCorrectionEvoPlot{fNo}, savename, extensions{k})
        end
    end
    % save (PeakPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\Peak_i.fig)
    for p = 1:numPeaks
        savename = [res_subfolder_path, '\', date, ...
            '_Peak_', num2str(p)];
        for k = 1:length(extensions)
            saveas(PeakPlot{p}, savename, extensions{k})
        end
    end

    % create subfolder for all Peak indicating plots
    res_subfolder_name = [date, '_Peak_Indication'];
    res_subfolder_path = [handles.datpath, '\', res_folder_name, '\', ...
        res_subfolder_name];
    if isfolder(res_subfolder_path) == 0
        mkdir(res_subfolder_path)
    end
    % save (FitSqPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\'TagName'.fig)
    for fNo = 1:fileNumber
        savename = [res_subfolder_path, '\', date, '_', ...
            TagName{fNo}];
        for k = 1:length(extensions)
            saveas(FitSqPlot{fNo}, ...
                strrep( strrep(savename,' ','_'), '.',''), extensions{k})
        end
    end
end
disp('S(q) analysis done');

%% ------------------------------------------------------------------------
% ------------------------------ COMPARE ----------------------------------
% -------------------------------------------------------------------------
function Button_Compare_Callback(hObject, eventdata, handles)
foldername = handles.datfname;
fileNumber = handles.fileNumber;
TagName    = evalin('base', 'TagName');
TagNameLegend = evalin('base', 'TagNameLegend');
TagTitle   = evalin('base', 'TagTitle');
mx   = handles.mx;
rmax = handles.rmax;
r   = evalin('base','r')';
Gr  = evalin('base','Gr')';
GrM = evalin('base','GrM')';

ref = 1;  % reference curve which is subtracted in second plot
PlotColors = [1;9;7.9];  % for cryst. PDF (Color(PlotColors))

%% import crystalline PDFs
[fnamePDF,pnamePDF] = uigetfile({'*.txt','TXT';'*.*','All files (*.*)';},...
    'Choose data file to compare with',...
    'Extended_eRDF_Anlyser',...
    'MultiSelect', 'on');

if iscell(fnamePDF) ~= 1
    if fnamePDF == 0
        % User clicked the Cancel button.
        return;
    end
    fnamePDF = char2cell(fnamePDF);
end

%% Plots

for PDFNo = 1:size(fnamePDF,2)

    PDFtxt  = importdata([pnamePDF,fnamePDF{PDFNo}], ',');
    contrib = PDFtxt.colheaders(1,2:end);
    Grt     = PDFtxt.data(:,2:end);
    rt      = PDFtxt.data(:,1);

    scal = 1;
    % rescale crystalline peak positions
    %     firstpeak_theo = 3; % 2.44;
    %     firstpeak_expe = 2.52;
    %     scal = round(firstpeak_expe/firstpeak_theo,3);
    rt = rt*scal;


    % find number of atom species
    Ntot       = size(Grt,2);
    NumSpecies = 0;
    n          = 1;
    while Ntot ~= 0
        Ntot       = Ntot - n;
        NumSpecies = NumSpecies + 1;
        n          = n + 1;
    end
    CombMixPDF = combnk(1:NumSpecies,2);

    % find atom species label
    species = cell(NumSpecies,1);
    for i = 1:NumSpecies
        ks = strfind(contrib(i), '-', 'ForceCellOutput', 1);
        ks = ks{1,1};
        species{i} = contrib{i}(1:ks-1);
    end


    % swap species entries if needed
    %     dummyspec  = species;
    %     species{1} = dummyspec{2};
    %     species{2} = dummyspec{1};
    %     dummyGrt = Grt;
    %     Grt(:,1)   = dummyGrt(:,2);
    %     Grt(:,2)   = dummyGrt(:,1);

    % ---------------------------------------------------------------------

    %% Plot experimental Gr and crystalline PDF
    ComparePlot = figure('Name', ['G(r): Compare to ', fnamePDF{PDFNo}], ...
        'units', 'normalized', 'outerposition', [0.1 0 0.7 1], ...
        'NumberTitle', 'off');

    if fileNumber == 1
        NumPlt = 1;
    else
        NumPlt = 2;
    end

    % set to 1 for Plot without DeltaG
    % ---|---
    % ---V---
    % NumPlt = 1;

    %tiledlayout(NumPlt,1)
    % sgtitle(TagTitle)
    ax1 = subplot(NumPlt,1,1);
    maxGr = 0;
    for fNo = 1:fileNumber
        if fNo == ref
            plot(ax1, r, GrM{fNo}, ...
                'DisplayName', [TagNameLegend{fNo}(:)', ' \textsf{($G_0$)}'],...
                'Color', Color(3+7.5/mx*(fNo-1)), ...
                'Linewidth', 1.0);
        else
            plot(ax1, r, GrM{fNo}, ...
                'DisplayName', TagNameLegend{fNo}(:), ...
                'Color', Color(3+7.5/mx*(fNo-1)), ...
                'Linewidth', 1.0);
        end
        hold on
        % find scaling
        rind = find(r == 1);
        if max(abs(GrM{fNo}(rind:end))) > maxGr
            maxGr = max(abs(GrM{fNo}(rind:end)));
        end
    end

    % plot reference line at Gr=zero
    xlims = get(handles.axes7, 'xlim');
    plot(ax1, [0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');

    % scale crystalline PDF
    Grt = Grt./max(sum(Grt,2)).*maxGr*1.05;

    %% Plot theoretical PDF

    % scale Icosahedron peaks to first maximum in Gr
    if string(fnamePDF{PDFNo}) == 'Icosahedron.txt'
        [Amx, Imx] = findpeaks(GrM{1,1});
        Amx2 = [];
        Imx2 = [];
        for i = 1:length(Imx)
            if Amx(i) > 0
                Amx2 = [Amx2, Amx(i)];
                Imx2 = [Imx2, Imx(i)];
            end
        end
        mImx = min(Imx2)./100;
    else
        mImx = 1.;
    end

    % plot
    b = bar(ax1, rt*mImx, Grt(:,:), 'stacked', 'barwidth', 3.0);

    % change plot colors and legend labels
    for k = 1:NumSpecies
        b(k).FaceColor   = element_color(species{k});
        b(k).EdgeColor   = 'k';
        b(k).LineWidth   = 0.2;
        % Plot notification if crystal structure has been scaled
        if scal ~= 1
            b(k).DisplayName = ['\textsf{', contrib{k}, ' (scaled by ', num2str(scal), ')}'];
        else
            b(k).DisplayName = ['\textsf{', contrib{k}, '}'];
        end

    end
    for nc = 1:size(CombMixPDF,1)
        k  = nc + NumSpecies;
        n1 = CombMixPDF(nc,1);
        n2 = CombMixPDF(nc,2);
        b(k).FaceColor   = element_color(species{n1});
        b(k).EdgeColor   = element_color(species{n2});
        b(k).LineWidth   = 1.3;
        b(k).LineStyle   = '--';

        b(k).DisplayName = ['\textsf{', contrib{k}, '}'];
    end
    xlabel(ax1, '$r$ [\AA]', 'Interpreter', 'LaTex');
    ylabel(ax1, '$G(r)$ [\AA$$^{-2}$$]', 'Interpreter', 'LaTex');
    %lgd = legend(ax1, 'show', 'Interpreter', 'Latex');
    xlim(ax1, [0 rmax]);

    % ---------------------------------------------------------------------

    %% Plot DeltaGr = Gr-Gr0 and crystalline PDF (if more than one curve)
    if NumPlt == 2
        maxDelta = 0;
        ax2 = subplot(2,1,2);
        for fNo = 1:fileNumber
            if fNo == ref
                plot(ax2, r, GrM{fNo}-GrM{ref}, ...
                    'DisplayName', [TagNameLegend{fNo}(:)', ' \textsf{($G_0$)}'], ...
                    'Color', Color(3+7.5/mx*(fNo-1)), ...
                    'HandleVisibility', 'off', ...
                    'Linewidth', 1.0);
            else
                plot(ax2, r, GrM{fNo}-GrM{ref}, ...
                    'DisplayName', TagNameLegend{fNo}(:), ...
                    'Color', Color(3+7.5/mx*(fNo-1)), ...
                    'HandleVisibility', 'off', ...
                    'Linewidth', 1.0);
                % non-cumulative difference
                %                 plot(ax2, r, GrM{fNo}-GrM{fNo-1}, ...
                %                     'DisplayName', TagNameLegend{fNo}(:),...
                %                     'Color', Color(3+7.5/mx*(fNo-1)));
            end
            hold on
            % find scaling
            rind = find(r == 1);
            if max(abs(GrM{fNo}(rind:end)-GrM{ref}(rind:end))) > maxDelta
                maxDelta = max(abs(GrM{fNo}(rind:end)-GrM{ref}(rind:end)));
            end
        end

        % scale crystalline PDF
        Grt = Grt./max(sum(Grt,2)).*maxDelta*1.05;

        % plot theoretical PDF
        % plot
        b = bar(ax2, rt*mImx, Grt(:,:), 'stacked', 'barwidth', 3.0);
        % change plot colors and legend labels
        for k = 1:NumSpecies
            b(k).FaceColor   = element_color(species{k});
            b(k).EdgeColor   = 'k';
            b(k).LineWidth   = 0.1;
            b(k).DisplayName = ['\textsf{', contrib{k}, '}'];
        end
        for nc = 1:size(CombMixPDF,1)
            k  = nc + NumSpecies;
            n1 = CombMixPDF(nc,1);
            n2 = CombMixPDF(nc,2);
            b(k).FaceColor   = element_color(species{n1});
            b(k).EdgeColor   = element_color(species{n2});
            b(k).LineWidth   = 1.3;
            b(k).LineStyle   = '--';

            b(k).DisplayName = ['\textsf{', contrib{k}, '}'];
        end
        xlabel(ax2, '$r$ [\AA]', 'Interpreter', 'LaTex');
        ylabel(ax2, '$\Delta G_i(r) = G_i(r) - G_0(r)$', ...
            'Interpreter', 'LaTex');
        xlim(ax2, [0 rmax]);
    else
        ax2 = ax1;
    end

    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
    lgd = legend(ax2, 'show', 'Interpreter', 'Latex');
    lgd.FontSize = 12;
    PrePos = lgd.Position;
    lgd.Position = [0.786944162378222, 0.804080481515227, PrePos(3), PrePos(4)];

    %% save figure if checkbox "save" is activated
    if handles.save == 1
        extensions = {'fig', 'svg', 'png'};
        rehash
        % create folder name (e.g. _02_03_04_05_06_Compare_to_GeTe_HT)
        res_folder_name = '_';
        for fNo = 1:fileNumber
            whitespace = find(isspace(TagName{fNo})==1, 1);
            res_folder_name = [res_folder_name, ...
                char(TagName{fNo}(1:whitespace-1)), ' '];
        end
        res_folder_name = [res_folder_name, 'Compare_to_', ...
            fnamePDF{PDFNo}(1:end-4)];
        res_folder_name = strrep(res_folder_name, ' ', '_');

        % delete latest figure if checkbox "overwrite" is activated
        % and if existing
        rehash
        if handles.overwrite == 1 && ...
                isfolder([handles.datpath, '\', res_folder_name]) == 1

            % find latest figure in folder
            FileInfo = dir([handles.datpath, '\', res_folder_name, '\*.fig']);
            TimeStamp = 0;
            delind = 0;
            for nfile = 1:size(FileInfo)
                FI = FileInfo(nfile,1);
                if FI.datenum > TimeStamp
                    TimeStamp = FI.datenum;
                    delind = nfile;
                end
            end

            % delete latest figure in folder
            FI = FileInfo(delind,1);
            for x = 1:length(extensions)
                delete([FI.folder, '\', FI.name(1:end-3), extensions{x}])
            end


            % create dir if not existing
        elseif isfolder([handles.datpath, '\', res_folder_name]) == 0
            mkdir([handles.datpath, '\', res_folder_name])
        end

        % get timestamp
        date = strrep(datestr(now), ':', '-');
        date = strrep(date, ' ', '_');

        % save the figure (02-Apr-2020_16-36-26.fig)
        savename = [handles.datpath, '\', res_folder_name, '\', ...
            date];
        for k = 1:length(extensions)
            saveas(ComparePlot, savename, extensions{k})
        end

    end

end

%% ------------------------------------------------------------------------
% ---------------------------- ANALYSE G(r) -------------------------------
% -------------------------------------------------------------------------
function Button_Analyse_Callback(hObject, eventdata, handles)
fileNumber = handles.fileNumber;
foldername = handles.datfname;
TagName    = evalin('base', 'TagName');
TagNameLegend = evalin('base', 'TagNameLegend');
TagTitle   = evalin('base', 'TagTitle');
Gr         = evalin('base', 'Gr');
GrM        = evalin('base', 'GrM');
CorrSteps  = evalin('base', 'CorrSteps');
numMean    = evalin('base', 'numMean');
rc         = evalin('base', 'rc');
CDegree    = evalin('base', 'CDegree');

rco          = evalin('base', 'rco');
rcoind       = evalin('base', 'rcoind');
MeanRangeInd = evalin('base', 'MeanRangeInd');

rmax = handles.rmax;
r    = 0.01:0.01:rmax;
mx   = handles.mx;

%[rminphys,~] = ginput(1);

rminphys = 2.2;
assignin('base', 'rminphys', rminphys);
% -------------------------------------------------------------------------

%% G(r) & Delta Gr = Gr-Gr0

ref = 1;  % index of reference curve

% Initialize Plot
GrPlot = figure('Name', 'G(r)', 'units', 'normalized', ...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');

% Plot G(r)
subplot(2,1,1);
grid on
hold on
box on
for fNo = 1:fileNumber
    if fNo == ref
        plot(r, GrM{fNo}, ...
            'DisplayName', [TagNameLegend{fNo}(:)', ' $$(G_0)$$'],...
            'Color', Color(3+7.5/mx*(fNo-1)), ...
            'Linewidth', 1.2);
    else
        plot(r, GrM{fNo}, ...
            'DisplayName', TagNameLegend{fNo}(:),...
            'Color', Color(3+7.5/mx*(fNo-1)), ...
            'Linewidth', 1.2);
    end
end
hold off
xlabel('$$r$$ [\AA]', 'Interpreter', 'Latex');
ylabel('$$G(r)$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
legend('show', 'Location', 'NorthEast', 'Interpreter', 'Latex');

% plot reference line at Gr=0
xlims = get(handles.axes7, 'xlim');
hold on
plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');
hold off

% Plot Delta Gr
subplot(2,1,2);
box on
% FFT = figure();
for fNo = 1:fileNumber

    DGrM = GrM{fNo}-GrM{ref};
    %     FFTDGr = fft(DGrM);
    %     Len = length(DGrM);
    %     P2 = abs(FFTDGr/Len);
    %     P1 = P2(1:Len/2+1);
    %     P1(2:end-1) = 2*P1(2:end-1);
    %     figure(FFT);
    %     hold on
    %     plot(P1);
    %     hold off
    %     figure(GrPlot);

    subplot(2,1,2);
    box on
    plot(r,DGrM, ...
        'DisplayName', TagNameLegend{fNo}(:),...
        'Color', Color(3+7.5/mx*(fNo-1)), ...
        'Linewidth', 1.3);
    hold on
end
xlabel('$$r$$ [\AA]', 'Interpreter', 'Latex');
ylabel('$$\Delta G_i(r) = G_i(r) - G_0(r)$$', 'Interpreter', 'Latex');


% Plot reference line at Delta Gr=0
xlims = get(handles.axes7, 'xlim');
ylims = get(gca, 'ylim');
ylim(1.1*ylims);
hold on
grid on
plot([1 xlims(2)], [0 0], 'k', 'HandleVisibility','off');
hold off

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

%% Fit Gauss to Peaks in Gr
numPeaks = 7; % peak = minimum OR maximum
rmin     = rminphys;
rmax     = max(r);
[~,rminind] = min(abs(r-rmin));
[~,rmaxind] = min(abs(r-rmax));

mu       = cell(1,fileNumber);
sigma    = cell(1,fileNumber);
scaling  = cell(1,fileNumber);
muE      = cell(1,fileNumber);
sigmaE   = cell(1,fileNumber);
scalingE = cell(1,fileNumber);

FitPlot    = cell(1,fileNumber);
GrEvoPlot  = cell(1,fileNumber);

% Computing fit
yl = [0 0];
for fNo = 1:fileNumber
    Gr0 = cellfun(@(x) x(rminind:rmaxind), Gr(:,fNo), 'UniformOutput', 0);
    r0  = r(rminind:rmaxind);

    %     figure('Name', 'S(q)-G(r) Evolution upon correction' ...
    %         + TagTitle + ' ' + TagName{fNo}, ...
    %         'units', 'normalized', ...
    %         'outerposition', [0.1 0 0.7 1], ...
    %         'NumberTitle', 'off');
    %     subplot(3,1,1)
    %     xlabel('$$r$$ [\AA]', 'Interpreter', 'Latex');
    %     ylabel('$$G(r)$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
    %     hold on
    %     box on
    %     Gr01 = cellfun(@(x) x, Gr(:,fNo), 'UniformOutput', 0);
    %         plot(r, Gr01{1}, 'Color', Color(3+7.5/(CorrSteps-1)*(1-1)), ...
    %             'DisplayName', '$G_\mathrm{u}(r)$ (uncorrected)', ...
    %             'Linewidth', 1.0);
    %     legend('Interpreter', 'Latex', 'Location', 'best')
    %
    %
    %
    %     subplot(3,1,2)
    %     xlabel('$$q$$ [\AA$$^{-1}$$]', 'Interpreter', 'Latex');
    %     ylabel('$$S(q)$$', 'Interpreter', 'Latex');
    %     hold on
    %     box on
    %     Sq01 = evalin('base', 'Sq01');
    %     for iC = 1:CorrSteps
    %         plot(q, Sq01{iC}, 'Color', Color(3+7.5/(CorrSteps-1)*(iC-1)), ...
    %             'DisplayName', ...
    %             sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA', ...
    %             rc(iC)), ...
    %             'Linewidth', 1.0);
    %     end
    %     % Plot reference line at Sq=0
    %     xlims = get(handles.axes6, 'xlim');
    %     hold on
    %     plot([0 xlims(2)], [1 1], 'k', 'HandleVisibility', 'off');
    %     hold off
    %     %legend('Interpreter', 'Latex', 'Location', 'southeast')
    %
    %
    %
    %     subplot(3,1,3)
    %     xlabel('$$r$$ [\AA]', 'Interpreter', 'Latex');
    %     ylabel('$$G(r)$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
    %     hold on
    %     box on
    %     Gr01 = cellfun(@(x) x, Gr(:,fNo), 'UniformOutput', 0);
    %     for iC = 1:CorrSteps
    %         plot(r, Gr01{iC}, 'Color', Color(3+7.5/(CorrSteps-1)*(iC-1)), ...
    %             'DisplayName', ...
    %             ... sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA$$\\,\\,$$($$C_\\mathrm{D}=%.2f$$)', ...
    %             ... rc(iC), CDegree(iC,fNo)), ...
    %             sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA', ...
    %             rc(iC)), ...
    %             'Linewidth', 1.0);
    %     end
    %     legend('Interpreter', 'Latex', 'Location', 'eastoutside', 'AutoUpdate', 'off')
    % %     chH = get(gca,'Children');  % bring first curve in foreground
    % %     set(gca,'Children',[chH(end);chH(1:end-1)])
    %     % same y limits
    %     cylim = get(gca, 'ylim');
    %     if cylim(1) < yl(1) && cylim(1) ~= 0
    %         yl(1) = cylim(1);
    %     end
    %     if cylim(2) > yl(2) && cylim(2) ~= 2
    %         yl(2) = cylim(2);
    %     end
    %     % plot reference line at Gr=0
    %     xlims = get(handles.axes7, 'xlim');
    %     hold on
    %     plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');
    %     hold off
    %
    %     % Set Fontsize
    %     set(findall(gcf,'-property','FontSize'), 'FontSize', 20)
    %
    %     [   mu{fNo},  sigma{fNo},  scaling{fNo}, ...
    %         muE{fNo}, sigmaE{fNo}, scalingE{fNo}, ...
    %         FitPlot{fNo}] = ...
    %         Gaussfit_PDF2( r0, Gr0, numPeaks, TagNameLegend{fNo}(:)' );
    %

    GrEvoPlot{fNo} = figure('Name', 'G(r) Evolution upon correction' ...
        + TagTitle + ' ' + TagName{fNo}, ...
        'units', 'normalized', ...
        'outerposition', [0.1 0 0.7 1], ...
        'NumberTitle', 'off');


    xlabel('$$r$$ [\AA]', 'Interpreter', 'Latex');
    ylabel('$$G(r,r_{\mathrm{co}})$$ [\AA$$^{-2}$$]', 'Interpreter', 'Latex');
    hold on
    box on
    Gr01 = cellfun(@(x) x, Gr(:,fNo), 'UniformOutput', 0);
    for iC = 1:CorrSteps
        plot(r, Gr01{iC}, 'Color', Color(3+7.5/(CorrSteps-1)*(iC-1)), ...
            'DisplayName', ...
            ... sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA$$\\,\\,$$($$C_\\mathrm{D}=%.2f$$)', ...
            ... rc(iC), CDegree(iC,fNo)), ...
            sprintf('$$r_\\mathrm{co} = %0.2f\\,$$\\AA', ...
            rc(iC)), ...
            'Linewidth', 1.0);
    end
    legend('Interpreter', 'Latex', 'Location', 'eastoutside', 'AutoUpdate', 'off')
    %     chH = get(gca,'Children');  % bring first curve in foreground
    %     set(gca,'Children',[chH(end);chH(1:end-1)])
    % same y limits
    cylim = get(gca, 'ylim');
    if cylim(1) < yl(1) && cylim(1) ~= 0
        yl(1) = cylim(1);
    end
    if cylim(2) > yl(2) && cylim(2) ~= 2
        yl(2) = cylim(2);
    end
    % plot reference line at Gr=0
    xlims = get(handles.axes7, 'xlim');
    hold on
    plot([0 xlims(2)], [0 0], 'k', 'HandleVisibility', 'off');
    hold off

    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

    [   mu{fNo},  sigma{fNo},  scaling{fNo}, ...
        muE{fNo}, sigmaE{fNo}, scalingE{fNo}, ...
        FitPlot{fNo}] = ...
        Gaussfit_PDF2( r0, Gr0, numPeaks, TagNameLegend{fNo}(:)' );

end
for fNo = 1:fileNumber
    figure(GrEvoPlot{fNo});
    ylim([yl(1) yl(2)]);
end

% -------------------------------------------------------------------------
% Sort out bad peaks / match peaks from different fits
% muS, sigmaS etc. only for the purpose of sorting (being deleted afterwards)
muS       = cell(length(mu),CorrSteps);
muES      = cell(length(mu),CorrSteps);
sigmaS    = cell(length(mu),CorrSteps);
sigmaES   = cell(length(mu),CorrSteps);
scalingS  = cell(length(mu),CorrSteps);
scalingES = cell(length(mu),CorrSteps);


% % first iteration: compare peaks of same correction degree iC
% for iC = MeanRangeInd
%
%     for fNo = 1:fileNumber
%         muS{fNo,iC}       = mu{fNo}(iC,:);
%         muES{fNo,iC}      = muE{fNo}(iC,:);
%         sigmaS{fNo,iC}    = sigma{fNo}(iC,:);
%         sigmaES{fNo,iC}   = sigmaE{fNo}(iC,:);
%         scalingS{fNo,iC}  = scaling{fNo}(iC,:);
%         scalingES{fNo,iC} = scalingE{fNo}(iC,:);
%     end
%
%     [~,refPeaksind] = max(cellfun('length', muS(:,iC)));
%     mrpInd          = min(refPeaksind);
%     refPeaks        = muS{mrpInd,iC};
%     fNo = 0;
%     while fNo < fileNumber
%         fNo = fNo + 1;
%         br  = 1;
%         i   = 0;
%         while (br == 1) && (i < length(muS{fNo,iC}))
%             i = i + 1;
%             if ~isnan(muS{fNo,iC}(i))
%                 [~,ind] = min(abs(refPeaks-muS{fNo,iC}(i)));
%                 if ind > i
%                     muS{mrpInd,iC} = ...
%                         [muS{mrpInd,iC}(1:i-1), muS{mrpInd,iC}(ind:end)];
%                     muES{mrpInd,iC} = ...
%                         [muES{mrpInd,iC}(1:i-1), muES{mrpInd,iC}(ind:end)];
%                     sigmaS{mrpInd,iC} = ...
%                         [sigmaS{mrpInd,iC}(1:i-1), sigmaS{mrpInd,iC}(ind:end)];
%                     sigmaES{mrpInd,iC} = ...
%                         [sigmaES{mrpInd,iC}(1:i-1), sigmaES{mrpInd,iC}(ind:end)];
%                     scalingS{mrpInd,iC} = ...
%                         [scalingS{mrpInd,iC}(1:i-1), scalingS{mrpInd,iC}(ind:end)];
%                     scalingES{mrpInd,iC} = ...
%                         [scalingES{mrpInd,iC}(1:i-1), scalingES{mrpInd,iC}(ind:end)];
%
%                     [~,refPeaksind] = max(cellfun('length',muS(:,iC)));
%                     mrpInd          = min(refPeaksind);
%                     refPeaks        = muS{mrpInd,iC};
%                     i   = length(muS{fNo,iC});
%                     fNo = 0;
%                     br  = 0;
%                 elseif ind < i
%                     if length(muS{fNo,iC}) > i
%                         muS{fNo,iC}       = [muS{fNo,iC}(1:ind),muS{fNo,iC}(i+1:end)];
%                         sigmaS{fNo,iC}    = [sigmaS{fNo,iC}(1:ind),sigmaS{fNo,iC}(i+1:end)];
%                         scalingS{fNo,iC}  = [scalingS{fNo,iC}(1:ind),scalingS{fNo,iC}(i+1:end)];
%                         muES{fNo,iC}      = [muES{fNo,iC}(1:ind),muES{fNo,iC}(i+1:end)];
%                         sigmaES{fNo,iC}   = [sigmaES{fNo,iC}(1:ind),sigmaES{fNo,iC}(i+1:end)];
%                         scalingES{fNo,iC} = [scalingES{fNo,iC}(1:ind),scalingES{fNo,iC}(i+1:end)];
%                     else
%                         muS{fNo,iC}       = muS{fNo,iC}(1:ind);
%                         sigmaS{fNo,iC}    = sigmaS{fNo,iC}(1:ind);
%                         scalingS{fNo,iC}  = scalingS{fNo,iC}(1:ind);
%                         muES{fNo,iC}      = muES{fNo,iC}(1:ind);
%                         sigmaES{fNo,iC}   = sigmaES{fNo,iC}(1:ind);
%                         scalingES{fNo,iC} = scalingES{fNo,iC}(1:ind);
%                     end
%                     i   = length(muS{fNo,iC});
%                     fNo = 0;
%                     br  = 0;
%                 elseif (i == length(muS{fNo,iC})) && (i < length(refPeaks))
%                     refPeaks = refPeaks(1:i);
%                     fNo = 0;
%                     br  = 0;
%                 end
%             else
%                 muS{fNo,iC}       = muS{fNo,iC}(1:i-1);
%                 sigmaS{fNo,iC}    = sigmaS{fNo,iC}(1:i-1);
%                 scalingS{fNo,iC}  = scalingS{fNo,iC}(1:i-1);
%                 muES{fNo,iC}      = muES{fNo,iC}(1:i-1);
%                 sigmaES{fNo,iC}   = sigmaES{fNo,iC}(1:i-1);
%                 scalingES{fNo,iC} = scalingES{fNo,iC}(1:i-1);
%             end
%         end
%     end
% end
% % second iteration: compare peaks of same fNo
% for fNo = 1:fileNumber
%
%     [~,refPeaksind] = max(cellfun('length', muS(fNo,:)));
%     mrpInd          = min(refPeaksind);
%     refPeaks        = muS{fNo,mrpInd};
%     iC = MeanRangeInd(1);
%     while iC < MeanRangeInd(end)
%         iC = iC + 1;
%         br = 1;
%         i  = 0;
%         while (br == 1) && (i < length(muS{fNo,iC}))
%             i = i + 1;
%             if ~isnan(muS{fNo,iC}(i))
%                 [~,ind] = min(abs(refPeaks-muS{fNo,iC}(i)));
%                 if ind > i
%                     muS{fNo,mrpInd} = ...
%                         [muS{fNo,mrpInd}(1:i-1), muS{fNo,mrpInd}(ind:end)];
%                     muES{fNo,mrpInd} = ...
%                         [muES{fNo,mrpInd}(1:i-1), muES{fNo,mrpInd}(ind:end)];
%                     sigmaS{fNo,mrpInd} = ...
%                         [sigmaS{fNo,mrpInd}(1:i-1), sigmaS{fNo,mrpInd}(ind:end)];
%                     sigmaES{fNo,mrpInd} = ...
%                         [sigmaES{fNo,mrpInd}(1:i-1), sigmaES{fNo,mrpInd}(ind:end)];
%                     scalingS{fNo,mrpInd} = ...
%                         [scalingS{fNo,mrpInd}(1:i-1), scalingS{fNo,mrpInd}(ind:end)];
%                     scalingES{fNo,mrpInd} = ...
%                         [scalingES{fNo,mrpInd}(1:i-1), scalingES{fNo,mrpInd}(ind:end)];
%
%                     [~,refPeaksind] = max(cellfun('length',muS(fNo,:)));
%                     mrpInd          = min(refPeaksind);
%                     refPeaks        = muS{fNo,mrpInd};
%                     i  = length(muS{fNo,iC});
%                     iC = 0;
%                     br = 0;
%                 elseif ind < i
%                     if length(muS{fNo,iC}) > i
%                         muS{fNo,iC}       = [muS{fNo,iC}(1:ind),muS{fNo,iC}(i+1:end)];
%                         sigmaS{fNo,iC}    = [sigmaS{fNo,iC}(1:ind),sigmaS{fNo,iC}(i+1:end)];
%                         scalingS{fNo,iC}  = [scalingS{fNo,iC}(1:ind),scalingS{fNo,iC}(i+1:end)];
%                         muES{fNo,iC}      = [muES{fNo,iC}(1:ind),muES{fNo,iC}(i+1:end)];
%                         sigmaES{fNo,iC}   = [sigmaES{fNo,iC}(1:ind),sigmaES{fNo,iC}(i+1:end)];
%                         scalingES{fNo,iC} = [scalingES{fNo,iC}(1:ind),scalingES{fNo,iC}(i+1:end)];
%                     else
%                         muS{fNo,iC}       = muS{fNo,iC}(1:ind);
%                         sigmaS{fNo,iC}    = sigmaS{fNo,iC}(1:ind);
%                         scalingS{fNo,iC}  = scalingS{fNo,iC}(1:ind);
%                         muES{fNo,iC}      = muES{fNo,iC}(1:ind);
%                         sigmaES{fNo,iC}   = sigmaES{fNo,iC}(1:ind);
%                         scalingES{fNo,iC} = scalingES{fNo,iC}(1:ind);
%                     end
%                     i  = length(muS{fNo,iC});
%                     iC = 0;
%                     br = 0;
%                 elseif (i == length(muS{fNo,iC})) && (i < length(refPeaks))
%                     refPeaks = refPeaks(1:i);
%                     iC = 0;
%                     br = 0;
%                 end
%             else
%                 muS{fNo,iC}       = muS{fNo,iC}(1:i-1);
%                 sigmaS{fNo,iC}    = sigmaS{fNo,iC}(1:i-1);
%                 scalingS{fNo,iC}  = scalingS{fNo,iC}(1:i-1);
%                 muES{fNo,iC}      = muES{fNo,iC}(1:i-1);
%                 sigmaES{fNo,iC}   = sigmaES{fNo,iC}(1:i-1);
%                 scalingES{fNo,iC} = scalingES{fNo,iC}(1:i-1);
%             end
%         end
%     end
% end
%
% % add results which do not contribute to averaging
% for iC = 1:(MeanRangeInd(1)-1)
%     for fNo = 1:fileNumber
%         muS{fNo,iC}       = mu{fNo}(iC,:);
%         muES{fNo,iC}      = muE{fNo}(iC,:);
%         sigmaS{fNo,iC}    = sigma{fNo}(iC,:);
%         sigmaES{fNo,iC}   = sigmaE{fNo}(iC,:);
%         scalingS{fNo,iC}  = scaling{fNo}(iC,:);
%         scalingES{fNo,iC} = scalingE{fNo}(iC,:);
%     end
% end
% for iC = (MeanRangeInd(end)+1):CorrSteps
%     for fNo = 1:fileNumber
%         muS{fNo,iC}       = mu{fNo}(iC,:);
%         muES{fNo,iC}      = muE{fNo}(iC,:);
%         sigmaS{fNo,iC}    = sigma{fNo}(iC,:);
%         sigmaES{fNo,iC}   = sigmaE{fNo}(iC,:);
%         scalingS{fNo,iC}  = scaling{fNo}(iC,:);
%         scalingES{fNo,iC} = scalingE{fNo}(iC,:);
%     end
% end
%
%
% % transfer matched peaks back to mu, sigma, etc.
% mu       = cell(1,fileNumber);
% sigma    = cell(1,fileNumber);
% scaling  = cell(1,fileNumber);
% muE      = cell(1,fileNumber);
% sigmaE   = cell(1,fileNumber);
% scalingE = cell(1,fileNumber);
% for fNo = 1:fileNumber
%     mu{fNo} = cell2mat(reshape(muS(fNo,:),[CorrSteps,1]));
%     muE{fNo} = cell2mat(reshape(muES(fNo,:),[CorrSteps,1]));
%     sigma{fNo} = cell2mat(reshape(sigmaS(fNo,:),[CorrSteps,1]));
%     sigmaE{fNo} = cell2mat(reshape(sigmaES(fNo,:),[CorrSteps,1]));
%     scaling{fNo} = cell2mat(reshape(scalingS(fNo,:),[CorrSteps,1]));
%     scalingE{fNo} = cell2mat(reshape(scalingES(fNo,:),[CorrSteps,1]));
% end
% clearvars muS muES sigmaS sigmaES scalingS scalingES
numPeaks = size(mu{1},2);

% Compute corrected averages of numMean correction ranges
muM       = cell(1,fileNumber);
sigmaM    = cell(1,fileNumber);
scalingM  = cell(1,fileNumber);
muEM      = cell(1,fileNumber);
sigmaEM   = cell(1,fileNumber);
scalingEM = cell(1,fileNumber);
for fNo = 1:fileNumber
    muM{fNo}       = mean(mu{fNo}(MeanRangeInd, :));
    muEM{fNo}      = std(mu{fNo}(MeanRangeInd, :));
    sigmaM{fNo}    = mean(sigma{fNo}(MeanRangeInd, :));
    sigmaEM{fNo}   = std(sigma{fNo}(MeanRangeInd, :));
    scalingM{fNo}  = mean(scaling{fNo}(MeanRangeInd, :));
    scalingEM{fNo} = std(scaling{fNo}(MeanRangeInd, :));
end

%% Define plot colors
cmu     = [0.00,0.45,0.74];   % green color for mu
csc     = Color(12.375);  % purpe color for scaling
csi     = Color(4.875);   % blue color for sigma

%% Plot correction evolution of Mu/Sigma/Scaling
numPlotPeaks = 3;  % number of peaks in the evolution plot (first 3 extrema)

PeakCorrectionEvoPlot = cell(fileNumber,1);
for fNo = 1:fileNumber
    PeakCorrectionEvoPlot{fNo} = figure('Name', ...
        'G(r) Mu/Sigma/Scaling Evolution upon correction: ' ...
        + TagTitle + ' ' + TagName{fNo}, ...
        'units', 'normalized', ...
        'outerposition', [0.1 0 0.7 1], ...
        'NumberTitle', 'off');
    sgtitle(TagNameLegend{fNo}, 'Interpreter', 'Latex');

    nmax = 0;
    nmin = 0;
    for p = 1:numPlotPeaks
        if scalingM{fNo}(end,p) >= 0
            nmax   = nmax + 1;
            extID  = 'Maximum';
            rsym   = 'r';
            sigsym = '\sigma';
            pID    = nmax;
        else
            nmin   = nmin + 1;
            extID  = 'Minimum';
            rsym   = '\tilde{r}';
            sigsym = '\tilde{\sigma}';
            pID    = nmin;
        end
        subplot(3,numPlotPeaks,p)
        hold on
        box on

        title([extID, ' ', num2str(pID)]);
        errorbar(mu{fNo}(:,p), muE{fNo}(:,p), 'Marker', 'x', ...
            'Linewidth', 1.0, 'DisplayName', '', ...
            'Color', cmu);
        plot_hline(muM{fNo}(p), muEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        xticklabels('')
        ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", rsym, pID), ...
            'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:3:end))
        set(gca,'XTick',1:3:CorrSteps)
        grid on

        subplot(3,numPlotPeaks,p+numPlotPeaks)
        hold on
        box on
        errorbar(sigma{fNo}(:,p), sigmaE{fNo}(:,p), 'Marker', 'x', ...
            'Linewidth', 1.0, 'DisplayName', '', ...
            'Color', csi)
        plot_hline(sigmaM{fNo}(p), sigmaEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        xticklabels({''})
        ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", sigsym, pID), ...
            'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:3:end))
        set(gca,'XTick',1:3:CorrSteps)
        grid on

        Ax = subplot(3,numPlotPeaks,p+2*numPlotPeaks);
        hold on
        box on
        errorbar(scaling{fNo}(:,p), scalingE{fNo}(:,p), 'Marker', 'x', ...
            'Linewidth', 1.0, 'DisplayName', '', ...
            'Color', csc)
        plot_hline(scalingM{fNo}(p), scalingEM{fNo}(p), ...
            rcoind-numMean-0.1, rcoind+numMean+0.1, ...
            'black')
        Ax.TickLabelInterpreter = 'none';
        ylabel(sprintf('$$G({%s}_{%i})$$', rsym, pID), ...
            'Interpreter', 'Latex')
        xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
        xlim([0.5 CorrSteps+0.5])
        set(gca,'XTickLabel',rc(1:3:end))
        set(gca,'XTick',1:3:CorrSteps)
        grid on
    end

    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
end

%% Plot correction evolution of r2/r1

r1r2CorrectionEvoPlot = figure('Name', ...
    'r2/r1 Evolution upon correction: ' ...
    + TagTitle + ' ' + TagName{fNo}, ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');

sgtitle('$r_2/r_1$ evolution upon correction', ...
    'Interpreter', 'Latex');
for fNo = 1:fileNumber
    r2r1   = mu{fNo}(:,3) ./ mu{fNo}(:,1);  % (:,2) is first minimum
    r2r1E  = sqrt( (muE{fNo}(:,3)./ mu{fNo}(:,3)).^2 ...
        + (muE{fNo}(:,1)./ mu{fNo}(:,1)).^2 ) .*abs(r2r1);
    r2r1M  = muM{fNo}(:,3) ./ muM{fNo}(:,1);
    r2r1M  = mean( r2r1(MeanRangeInd) );
    r2r1EM = sqrt( (muEM{fNo}(:,3)./ muM{fNo}(:,3)).^2 ...
        + (muEM{fNo}(:,1)./ muM{fNo}(:,1)).^2 ) .*abs(r2r1M);
    r2r1EM =  std( r2r1(MeanRangeInd) );

    %Ax = subplot(fileNumber,1,fNo);
    Ax = subplot(1,1,1);
    hold on
    box on
    errorbar(r2r1,r2r1E, '-', 'Marker', 'x', ...
        'DisplayName', TagNameLegend{fNo}(:),...
        'Color', Color(3+7.5/mx*(fNo-1)),...
        'Linewidth', 1.0);
    plot_hline(r2r1M, r2r1EM, ...
        rcoind-numMean-0.1, rcoind+numMean+0.1, ...
        Color(3+7.5/mx*(fNo-1)))

    Ax.TickLabelInterpreter = 'none';
    ylabel('$r_2/r_1$', 'Interpreter', 'Latex')
    xlabel('$$r_\mathrm{co}$$ [\AA]', 'Interpreter', 'Latex')
    xlim([0.5 CorrSteps+0.5])
    set(gca,'XTickLabel',rc)
    set(gca,'XTick',1:CorrSteps)

    grid on
    legend('show', 'Location', 'SouthEast', 'Interpreter', 'Latex');

    % same y limits
    if fNo == 1
        yl = [min(r2r1) max(r2r1)];
    end
    cylim = get(gca, 'ylim');
    if cylim(1) < yl(1)
        yl(1) = cylim(1);
    end
    if cylim(2) > yl(2)
        yl(2) = cylim(2);
    end
end

% same scale for y-axes
% for fNo = 1:fileNumber
%     subplot(fileNumber,1,fNo);
%     ylim([yl(1) yl(2)]);
% end

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
% -------------------------------------------------------------------------

%% Error propagation
for fNo = 1:fileNumber
    for p = 1:numPeaks
        muCorrErr    = muEM{fNo}(p);
        muFitErr     = max(muE{fNo}(MeanRangeInd,p));
        muEM{fNo}(p) = sqrt( muCorrErr^2 + muFitErr^2 );
        % estimating the relative error on mu to 0.05% (from results of mu from different SAEDs of same sample/annealing state)
        muEM{fNo}(p) = sqrt( muEM{fNo}(p)^2 + (0.0005.*muM{fNo}(p))^2 );

        sigmaCorrErr    = sigmaEM{fNo}(p);
        sigmaFitErr     = max(sigmaE{fNo}(MeanRangeInd,p));
        sigmaEM{fNo}(p) = sqrt( sigmaCorrErr^2 + sigmaFitErr^2 );

        scalingCorrErr    = scalingEM{fNo}(p);
        scalingFitErr     = max(scalingE{fNo}(MeanRangeInd,p));
        scalingEM{fNo}(p) = sqrt( scalingCorrErr^2 + scalingFitErr^2 );
    end
end

%% Re-organize results for plotting
muM        = cell2mat(muM');
sigmaM     = cell2mat(sigmaM');
scalingM   = cell2mat(scalingM');
muEM       = cell2mat(muEM');
sigmaEM    = cell2mat(sigmaEM');
scalingEM  = cell2mat(scalingEM');

%% Plotting mu, (sigma), scaling
PeakPlot  = cell(fileNumber,1);
nmax = 0;
nmin = 0;

for p = 1:numPeaks
    if scalingM(1,p) >= 0
        nmax   = nmax + 1;
        extID  = 'Maximum';
        rsym   = 'r';
        sigsym = '\sigma';
        pID    = nmax;
    else
        nmin   = nmin + 1;
        extID  = 'Minimum';
        rsym   = '\tilde{r}';
        sigsym = '\tilde{\sigma}';
        pID    = nmin;
    end

    PeakPlot{p} = figure('Name', 'G(r) Mu/(Sigma/)Scaling: ' ...
        + TagTitle + [': G(r) Local extremum ', ...
        int2str(p), ' (', extID, ')'], ...
        'units','normalized', ...
        'outerposition',[0.1 0 0.7 1], ...
        'NumberTitle', 'off');

    % mu
    subplot(2,1,1)
    sgtitle(['$$G(r)$$ ', extID, ' ', int2str(pID)]);
    errorbar(muM(:,p), muEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'Linewidth', 1.0, ...
        'Color', cmu)
    xticklabels('')
    ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", rsym, pID), ...
        'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(gca,'XTick',1:numel(TagName))
    grid on

    % sigma
    %     subplot(3,1,2)
    %     errorbar(sigmaM(:,p), sigmaEM(:,p), 'Marker', 'x', ...
    %         'DisplayName', '', 'Linewidth', 1.0, ...
    %         'Color', csi)
    %     xticklabels({''})
    %     ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", sigsym, pID), ...
    %         'Interpreter', 'Latex')
    %     xlim([0.5 fileNumber+0.5])
    %     set(gca,'XTick',1:numel(TagName))
    %     grid on

    % scaling
    Ax = subplot(2,1,2);
    errorbar(scalingM(:,p), scalingEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'Linewidth', 1.0, ...
        'Color', csc)
    Ax.TickLabelInterpreter = 'none';
    ylabel(sprintf('$$G({%s}_{%i})$$', rsym, pID), ...
        'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'XTickLabel',TagNameLegend)
    set(gca,'XTick',1:numel(TagNameLegend))
    grid on


    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
end

%% Plotting all Mu
AllMuPlot = figure('Name', 'G(r) Mu', ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');
tx = tiledlayout(ceil(numPeaks/2),2);
tx.TileSpacing = 'compact';
tx.Padding     = 'compact';

nmax = 0;
nmin = 0;
for p = 1:numPeaks
    if scalingM(1,p) >= 0
        nmax   = nmax + 1;
        extID  = 'Maxima';
        rsym   = 'r';
        pID    = nmax;
    else
        nmin   = nmin + 1;
        extID  = 'Minima';
        rsym   = '\tilde{r}';
        pID    = nmin;
    end

    nexttile
    errorbar(muM(:,p), muEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'Linewidth', 1.0, ...
        'Color', cmu);

    if p == 1 || p == 2
        title(extID);
    end

    xticklabels('')
    ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", rsym, pID), ...
        'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(gca,'XTick',1:numel(TagNameLegend))
    if p == numPeaks || p == numPeaks-1
        Ax = gca;
        set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
        set(gca,'XTickLabel',TagNameLegend)
        xtickangle(11)
    end
    box on
end

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

%% Plotting all Scaling

AllScalingPlot = figure('Name', 'G(r) Scaling', ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');
tx = tiledlayout(ceil(numPeaks/2),2);
tx.TileSpacing = 'compact';
tx.Padding     = 'compact';

nmax = 0;
nmin = 0;
for p = 1:numPeaks
    if scalingM(1,p) >= 0
        nmax   = nmax + 1;
        extID  = 'Maxima';
        rsym   = 'r';
        pID    = nmax;
    else
        nmin   = nmin + 1;
        extID  = 'Minima';
        rsym   = '\tilde{r}';
        pID    = nmin;
    end

    nexttile
    errorbar(scalingM(:,p), scalingEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'Linewidth', 1.0, ...
        'Color', csc);

    if p == 1 || p == 2
        title(extID);
    end

    xticklabels('')
    ylabel(sprintf("$$G({%s}_{%i})$$", rsym, pID), ...
        'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(gca,'XTick',1:numel(TagNameLegend))
    if p == numPeaks || p == numPeaks-1
        Ax = gca;
        set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
        set(gca,'XTickLabel',TagNameLegend)
        xtickangle(11)
    end
    box on
end

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

%% Plotting all Sigma

AllSigmaPlot = figure('Name', 'G(r) Sigma', ...
    'units','normalized', ...
    'outerposition',[0.1 0 0.7 1], ...
    'NumberTitle', 'off');
tx = tiledlayout(ceil(numPeaks/2),2);
tx.TileSpacing = 'compact';
tx.Padding     = 'compact';

nmax = 0;
nmin = 0;
for p = 1:numPeaks
    if scalingM(1,p) >= 0
        nmax   = nmax + 1;
        extID  = 'Maxima';
        sigsym = '\sigma';
        pID    = nmax;
    else
        nmin   = nmin + 1;
        extID  = 'Minima';
        sigsym = '\tilde{\sigma}';
        pID    = nmin;
    end

    nexttile
    errorbar(sigmaM(:,p), sigmaEM(:,p), 'Marker', 'x', ...
        'DisplayName', '', 'Linewidth', 1.0, ...
        'Color', csi);

    if p == 1 || p == 2
        title(extID);
    end

    xticklabels('')
    ylabel(sprintf("$${%s}_{%i}$$ [\\AA]", sigsym, pID), ...
        'Interpreter', 'Latex')
    xlim([0.5 fileNumber+0.5])
    set(gca,'XTick',1:numel(TagNameLegend))
    if p == numPeaks || p == numPeaks-1
        Ax = gca;
        set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
        set(gca,'XTickLabel',TagNameLegend)
        xtickangle(11)
    end
    box on
end

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

%% Plotting r1, r2, r2/r1

% Find first 2 maxima
maxima = (scalingM > 0);  % true for all positive peaks (maxima)
maxima = maxima(1,:);
p = cell(2,1);
c = 1;
i = 1;
while (c <= 2) && (i <= size(maxima,2))  % c counts the first 2 maxima
    if maxima(i) == 1
        p{c,1} = i;
        c = c + 1;
    end
    i = i + 1;
end
p1 = p{1,1};
p2 = p{2,1};

% Plot
r1r2Plot = figure('Name', 'G(r) Peak Shift ' + TagTitle,...
    'units', 'normalized', ...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');

subplot(3,1,1)
%sgtitle(TagTitle + ': $r_2/r_1$', 'Interpreter', 'Latex');
errorbar(muM(:,p1), muEM(:,p1), 'Marker', 'x', ...
    'DisplayName', '', 'Linewidth', 1.0)
xticklabels({''})
ylabel('$$r_1$$ [\AA]', 'Interpreter', 'Latex')
set(gca,'XTick',1:numel(TagName))
xlim([0.5 fileNumber+0.5])
grid on

subplot(3,1,2)
errorbar(muM(:,p2), muEM(:,p2), 'Marker', 'x', ...
    'DisplayName', '', 'Linewidth', 1.0)
xticklabels({''})
ylabel('$r_2$ [\AA]', 'Interpreter', 'Latex')
set(gca,'XTick',1:numel(TagName))
xlim([0.5 fileNumber+0.5])
grid on


% r2/r1
r2r1  = muM(:,p2) ./ muM(:,p1);
r2r1E = sqrt( (muEM(:,p2)./ muM(:,p2)).^2 ...
    + (muEM(:,p1)./ muM(:,p1)).^2 ).*abs(r2r1);

Ax = subplot(3,1,3);
errorbar(r2r1, r2r1E, 'Marker', 'x', ...
    'DisplayName', '', 'Linewidth', 1.0)
Ax.TickLabelInterpreter = 'none';
ylabel('$r_2/r_1$', 'Interpreter', 'Latex')
set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
set(gca,'XTickLabel',TagNameLegend)
set(gca,'XTick',1:numel(TagNameLegend))
xtickangle(8)
xlim([0.5 fileNumber+0.5])
grid on

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
% -------------------------------------------------------------------------

%% Plot average bond angle

alpha  = 2*asind(r2r1/2);
alphaE = 1./sqrt(1-(r2r1./2).^2).*r2r1E*180/pi;

% Plot
AlphaPlot = figure('Name', 'Bond Angle Shift ' + TagTitle, ...
    'units', 'normalized', ...
    'outerposition', [0.1 0 0.5 1], ...
    'NumberTitle', 'off');

Ax = subplot(1,1,1);
errorbar(alpha, alphaE, 'Marker', 'x', ...
    'HandleVisibility', 'off', 'Linewidth', 1.0)
Ax.TickLabelInterpreter = 'none';
ylabel('\textsf{Average Bond Angle } $\alpha$', 'Interpreter', 'Latex')
set(Ax.XAxis,'TickLabelInterpreter', 'Latex');
set(gca,'XTickLabel',TagNameLegend)
set(gca,'XTick',1:numel(TagNameLegend))
xtickangle(8)
xlim([0.5 fileNumber+0.5])
grid on

% Plot reference line at alpha=90, 96, 109
xlims = get(gca, 'xlim');
hold on
grid on

xs  = linspace(xlims(1), xlims(2), 15);
ys1 = linspace(109.5, 109.5, 15);  % tetrahedral
ys2GeTe = linspace(95.82, 95.82, 15);  % GeTe PD (most prominent, r_2(short)/r_1(short))
ys2GeSe = linspace(95.2749, 95.2749, 15);  % GeSe PD
ys2Sb2Te = linspace(96.1, 96.1, 15);  % Sb2Te PD
ys3 = linspace(90, 90, 15);  % cubic

ys2 = ys2GeTe;

plot(xs, ys1, '--k^', ...
    'DisplayName', 'Tetrahedral');
plot(xs, ys2, '--kp', ...
    'DisplayName', 'Peierls distorted cubic (GeTe)');
% plot(xs, ys2, '--kp', ...
%     'DisplayName', 'Peierls distorted cubic (GeSe)');
plot(xs, ys3, '--ks', ...
    'DisplayName', 'Cubic');

hold off
legend('Location', 'best')
ylim([89.5 110.5])


% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

%% Create 3D plot if fileNumber > 1
if fileNumber > 1
    % ---------------------------------------------------------------------
    % Peak Positions (PP)
    MeshPlotPP = figure('Name', ...
        'G(r) 3D Plot: Peak Positions ' + TagTitle, ...
        'units', 'normalized', ...
        'outerposition', [0 0 1 1], ...
        'NumberTitle', 'off');
    p = mean(muM,1);
    refmu = muM(1,:);
    fNo = 1:(numel(TagName));

    relShiftPP = (muM-refmu)./refmu;

    t1             = tiledlayout(1,2);
    t1.TileSpacing = 'compact';
    t1.Padding     = 'compact';



    nexttile();
    mesh(p, fNo, relShiftPP, 'FaceAlpha', '0.4', ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on

    % plot contour lines
    contour3(p, fNo, relShiftPP, 4, '--k');
    contour3(p, fNo, relShiftPP, [0 0], '-k');

    % set labels etc.
    %title(TagTitle);
    colorbar;
    xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    xlim([min(p) max(p)]);
    xticks(round(p,3));
    xtickangle(-30)
    ylim([1 numel(TagName)]);
    Ax = gca;
    set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'YTickLabel',TagNameLegend)
    set(gca,'YTick',1:numel(TagNameLegend))



    ytickangle(30)
    zlabel('$$\frac{r_i-r_0}{r_0}$$', 'Interpreter', 'Latex', ...
        'Fontsize', 14);

    % plot reference Gr
    refG     = GrM{1}(rminind:rmaxind);
    refr     = r(rminind:rmaxind);
    maxG     = max(max( abs(refG) ));
    maxrel   = max(max( abs(relShiftPP) ));
    G_scaled = 0.75*refG./maxG*maxrel;
    yPlane   = 0*refr+1;
    plot3(refr, yPlane, G_scaled, 'k');


    % -----------------------------
    nexttile();
    mesh(p, fNo, relShiftPP, 'FaceAlpha', '0.4', ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on
    box on

    % plot contour lines
    contour3(p, fNo, relShiftPP, 4, '--k');
    contour3(p, fNo, relShiftPP, [0 0], '-k');

    % set labels etc.
    %title(TagTitle);
    xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    xlim([min(p) max(p)]);
    xticks(round(p,3));
    xtickangle(-30)
    ylim([0 numel(TagName)]);
    Ax = gca;
    set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'YTickLabel',TagNameLegend)
    set(gca,'YTick',1:numel(TagNameLegend))
    ytickangle(35)


    % plot flat reference Gr
    Gspan         = max(refG)+abs(min(refG));
    G_scaled_flat = refG+abs(min(refG));
    G_scaled_flat = G_scaled_flat./Gspan;
    zPlane        = 0*refr;
    text(p(1)+0.7*(p(2)-p(1)), 0.6, '$G(r)$ [\AA$^{-2}$]')
    plot3(refr, G_scaled_flat, zPlane, 'k');
    view([0 0 1]);


    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
    % ---------------------------------------------------------------------
    % Peak Heigths (PH)
    MeshPlotPH = figure('Name', ...
        'G(r) 3D Plot: Peak Heights ' + TagTitle, ...
        'units', 'normalized', ...
        'outerposition', [0 0 1 1], ...
        'NumberTitle', 'off');
    p = mean(muM,1);
    refsca = scalingM(1,:);
    fNo = 1:(numel(TagName));

    relShiftPH = (scalingM-refsca)./refsca;

    t2             = tiledlayout(1,2);
    t2.TileSpacing = 'compact';
    t2.Padding     = 'compact';


    nexttile();
    mesh(p, fNo, relShiftPH, 'FaceAlpha', '0.4', ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on

    % plot contour lines
    contour3(p, fNo, relShiftPH, 3, '--k');
    contour3(p, fNo, relShiftPH, [0 0], '-k');

    % set labels etc.
    %title(TagTitle);
    colorbar;
    xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    xlim([min(p) max(p)]);
    xticks(round(p,3));
    xtickangle(-30)
    ylim([1 numel(TagName)]);
    Ax = gca;
    set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'YTickLabel',TagNameLegend)
    set(gca,'YTick',1:numel(TagNameLegend))
    ytickangle(30)
    zlabel('$$\frac{G(r_i)-G(r_0)}{G(r_0)}$$', 'Interpreter', 'Latex', ...
        'Fontsize', 14);

    % plot reference Gr
    refG     = GrM{1}(rminind:rmaxind);
    refr     = r(rminind:rmaxind);
    maxG     = max(max( abs(refG) ));
    maxrel   = max(max( abs(relShiftPH) ));
    G_scaled = 0.75*refG./maxG*maxrel;
    yPlane   = 0*refr+1;
    plot3(refr, yPlane, G_scaled, 'k');


    % -----------------------------
    nexttile();
    mesh(p, fNo, relShiftPH, 'FaceAlpha', '0.4', ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on
    box on

    % plot contour lines
    contour3(p, fNo, relShiftPH, 3, '--k');
    contour3(p, fNo, relShiftPH, [0 0], '-k');

    % set labels etc.
    %title(TagTitle);
    xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    xlim([min(p) max(p)]);
    xticks(round(p,3));
    xtickangle(-30)
    ylim([0 numel(TagName)]);
    Ax = gca;
    set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    set(gca,'YTickLabel',TagNameLegend)
    set(gca,'YTick',1:numel(TagNameLegend))
    ytickangle(35)

    % plot flat reference Gr
    Gspan         = max(refG)+abs(min(refG));
    G_scaled_flat = refG+abs(min(refG));
    G_scaled_flat = G_scaled_flat./Gspan;
    zPlane        = 0*refr;
    text(p(1)+0.7*(p(2)-p(1)), 0.6, '$G(r)$ [\AA$^{-2}$]')
    plot3(refr, G_scaled_flat, zPlane, 'k');
    view([0 0 1]);


    % Set Fontsize
    set(findall(gcf,'-property','FontSize'), 'FontSize', 14)

    %     % ---------------------------------------------------------------------
    %     % Overall Shift (OS)
    %     MeshPlotOS = figure('Name', ...
    %         'G(r) 3D Plot: Overall Shift ' + TagTitle, ...
    %         'units', 'normalized', ...
    %         'outerposition', [0 0 1 1], ...
    %         'NumberTitle', 'off');
    %     p = mean(muM,1);
    %     fNo = 1:(numel(TagName));
    %
    %     NormPHPP = max(max( abs(relShiftPP) )) ./ max(max( abs(relShiftPH) ));
    %     % relShiftOS = sqrt( (relShiftPP.*NormPHPP).^2 + relShiftPH.^2 );
    %
    %     relShiftOS = sqrt( (relShiftPP./max(max( abs(relShiftPP) ))).^2 ...
    %                     +  (relShiftPH./max(max( abs(relShiftPH) ))).^2 );
    %
    %     t3             = tiledlayout(1,2);
    %     t3.TileSpacing = 'compact';
    %     t3.Padding     = 'compact';
    %
    %     nexttile();
    %     mesh(p, fNo, relShiftOS, 'FaceAlpha', '0.4', ...
    %         'FaceColor', 'interp', 'EdgeColor', 'interp');
    %     hold on
    %
    %     % plot contour lines
    %     contour3(p, fNo, relShiftOS, 3, '--k');
    %     contour3(p, fNo, relShiftOS, [0 0], '-k');
    %
    %     % set labels etc.
    %     %title(TagTitle);
    %     colorbar;
    %     xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    %     xlim([min(p) max(p)]);
    %     xticks(round(p,3));
    %     xtickangle(-30)
    %     ylim([1 numel(TagName)]);
    %     Ax = gca;
    %     set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    %     set(gca,'YTickLabel',TagNameLegend)
    %     set(gca,'YTick',1:numel(TagNameLegend))
    %     ytickangle(30)
    %     zlabel('Relative Overall-Shift', 'Interpreter', 'Latex')
    %
    %     % plot reference Gr
    %     refG     = GrM{1}(rminind:rmaxind);
    %     refr     = r(rminind:rmaxind);
    %     maxG     = max(max( abs(refG) ));
    %     maxrel   = max(max( abs(relShiftOS) ));
    %     G_scaled = 0.75*refG./maxG*maxrel;
    %     yPlane   = 0*refr+1;
    %     plot3(refr, yPlane, G_scaled, 'k');
    %
    %     nexttile();
    %     mesh(p, fNo, relShiftOS, 'FaceAlpha', '0.4', ...
    %         'FaceColor', 'interp', 'EdgeColor', 'interp');
    %     hold on
    %     box on
    %
    %     % plot contour lines
    %     contour3(p, fNo, relShiftOS, 3, '--k');
    %     contour3(p, fNo, relShiftOS, [0 0], '-k');
    %
    %     % set labels etc.
    %     %title(TagTitle);
    %     xlabel('Average Extremum Position [\AA]', 'Interpreter', 'Latex')
    %     xlim([min(p) max(p)]);
    %     xticks(round(p,3));
    %     xtickangle(-30)
    %     ylim([0 numel(TagName)]);
    %     Ax = gca;
    %     set(Ax.YAxis,'TickLabelInterpreter', 'Latex');
    %     set(gca,'YTickLabel',TagNameLegend)
    %     set(gca,'YTick',1:numel(TagNameLegend))
    %     ytickangle(35)
    %
    %     % plot flat reference Gr
    %     Gspan         = max(refG)+abs(min(refG));
    %     G_scaled_flat = refG+abs(min(refG));
    %     G_scaled_flat = G_scaled_flat./Gspan;
    %     zPlane        = 0*refr;
    %     text(p(1)+0.7*(p(2)-p(1)), 0.6, '$G(r)$ [\AA$^{-2}$]')
    %     plot3(refr, G_scaled_flat, zPlane, 'k');
    %     view([0 0 1]);
    %
    %
    %     % Set Fontsize
    %     set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
end

% -------------------------------------------------------------------------

%% save figures if checkbox "save" is activated
if handles.save == 1
    disp('saving...')
    extensions = {'fig', 'svg', 'png'};
    rehash
    % create folder name (e.g. _02_03_04_05_06_Analyse_G(r))
    res_folder_name = '_';
    for fNo = 1:fileNumber
        whitespace = find(isspace(TagName{fNo})==1, 1);
        res_folder_name = [res_folder_name, ...
            char(TagName{fNo}(1:whitespace-1)), ' '];
    end
    res_folder_name = [res_folder_name, 'Analyse_G(r)'];
    res_folder_name = strrep(res_folder_name, ' ', '_');

    % delete latest data if checkbox "overwrite" is activated
    % and if existing
    rehash
    if handles.overwrite == 1 && ...
            isfolder([handles.datpath, '\', res_folder_name]) == 1

        % find latest figure in folder
        FileInfo = dir([handles.datpath, '\', res_folder_name, '\*.fig']);
        TimeStamp = 0;
        delind = 0;
        for nfile = 1:size(FileInfo)
            FI = FileInfo(nfile,1);
            if FI.datenum > TimeStamp
                TimeStamp = FI.datenum;
                delind = nfile;
            end
        end

        % find all folders and files with same date
        FI = FileInfo(delind,1);
        sim_name = FI.name(1:20);
        for k = 1:length(extensions)
            delete([FI.folder, '\', sim_name, '*.', extensions{k}]);
        end
        delfolders = dir([FI.folder, '\', sim_name, '*']);
        for i = 1:size(delfolders)
            deldata = delfolders(i);
            delpath = [deldata.folder, '\', deldata.name];
            rmdir(delpath, 's');
        end


        % create dir if not existing
    elseif isfolder([handles.datpath, '\', res_folder_name]) == 0
        mkdir([handles.datpath, '\', res_folder_name])
    end

    % get timestamp
    date = strrep(datestr(now), ':', '-');
    date = strrep(date, ' ', '_');


    % save (GrPlot) figure (02-Apr-2020_16-36-26_Gr.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_Gr'];
    for k = 1:length(extensions)
        saveas(GrPlot, savename, extensions{k})
    end
    % save (GrEvoPlot) figure (02-Apr-2020_16-36-26_Gr_Correction_Evolution.fig)
    for fNo = 1:fileNumber
        savename = [handles.datpath, '\', res_folder_name, '\', ...
            date, '_Gr_Correction_Evolution_', strrep(TagName{fNo}, '.','')];
        for k = 1:length(extensions)
            saveas(GrEvoPlot{fNo}, savename, extensions{k})
        end
    end


    % save r1r2CorrectionEvoPlot (02-Apr-2020_16-36-26_r1r2_Correction_Evolution.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_r1r2_Correction_Evolution'];
    for k = 1:length(extensions)
        saveas(r1r2CorrectionEvoPlot, savename, extensions{k})
    end
    % save r1r2Plot (02-Apr-2020_16-36-26_r1r2.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_r1r2'];
    for k = 1:length(extensions)
        saveas(r1r2Plot, savename, extensions{k})
    end
    % save AlphaPlot (02-Apr-2020_16-36-26_AlphaPlot.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_AlphaPlot'];
    for k = 1:length(extensions)
        saveas(AlphaPlot, savename, extensions{k})
    end

    % save AllMuPlot (02-Apr-2020_16-36-26_All_Mu.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_All_Mu'];
    for k = 1:length(extensions)
        saveas(AllMuPlot, savename, extensions{k})
    end
    % save AllScalingPlot (02-Apr-2020_16-36-26_All_Scaling.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_All_Scaling'];
    for k = 1:length(extensions)
        saveas(AllScalingPlot, savename, extensions{k})
    end
    % save AllSigmaPlot (02-Apr-2020_16-36-26_All_Sigma.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_All_Sigma'];
    for k = 1:length(extensions)
        saveas(AllSigmaPlot, savename, extensions{k})
    end

    % save MeshPlotPP (02-Apr-2020_16-36-26_MeshPlot_Peak_Position_Shifts.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_MeshPlot_Peak_Position_Shifts'];
    for k = 1:length(extensions)
        saveas(MeshPlotPP, savename, extensions{k})
    end
    % save MeshPlotPH (02-Apr-2020_16-36-26_MeshPlot_Peak_Height_Shifts.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date, '_MeshPlot_Peak_Height_Shifts'];
    for k = 1:length(extensions)
        saveas(MeshPlotPH, savename, extensions{k})
    end
%     % save MeshPlotOS (02-Apr-2020_16-36-26_MeshPlot_Overall_Shifts.fig)
%     savename = [handles.datpath, '\', res_folder_name, '\', ...
%         date, '_MeshPlot_Overall_Shifts'];
%     for k = 1:length(extensions)
%         saveas(MeshPlotOS, savename, extensions{k})
%     end


    % create subfolder for all Peakfit results
    res_subfolder_name = [date, '_Peakfit_Results'];
    res_subfolder_path = [handles.datpath, '\', res_folder_name, '\', ...
        res_subfolder_name];
    if isfolder(res_subfolder_path) == 0
        mkdir(res_subfolder_path)
    end
    % save (PeakCorrectionEvoPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\Evolution_Extremum_1+2+3.fig)
    for fNo = 1:fileNumber
        savename = [res_subfolder_path, '\', date, '_', ...
            strrep( strrep(TagName{fNo},' ','_'), '.',''), ...
            '_Evolution_Extremum_1+2+3'];
        for k = 1:length(extensions)
            saveas(PeakCorrectionEvoPlot{fNo}, savename, extensions{k})
        end
    end
    % save (PeakPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\Peak_i.fig)
    for p = 1:numPeaks
        savename = [res_subfolder_path, '\', date, ...
            '_Peak_', num2str(p)];
        for k = 1:length(extensions)
            saveas(PeakPlot{p}, savename, extensions{k})
        end
    end


    % create subfolder for all Peak indicating plots
    res_subfolder_name = [date, '_Peak_Indication'];
    res_subfolder_path = [handles.datpath, '\', res_folder_name, '\', ...
        res_subfolder_name];
    if isfolder(res_subfolder_path) == 0
        mkdir(res_subfolder_path)
    end
    % save (FitPlot) figures (02-Apr-2020_16-36-26_Peakfit_Results\'TagName'.fig)
    for fNo = 1:fileNumber
        savename = [res_subfolder_path, '\', date, '_', ...
            TagName{fNo}];
        for k = 1:length(extensions)
            saveas(FitPlot{fNo}, strrep(strrep(savename, ' ','_'), '.', ','), extensions{k})
        end
    end
end
disp('G(r) analysis done');



%% ------------------------------------------------------------------------
% ---------------------------- EXPORT (NEW) -------------------------------
% -------------------------------------------------------------------------
function Button_Export_Fig_Callback(hObject, eventdata, handles)
TagTitle = evalin('base', 'TagTitle');
TagName = evalin('base', 'TagName');
fileNumber = handles.fileNumber;

%% Prepare export data for Sq & Gr
SqM     = evalin('base', 'SqM');
GrM     = evalin('base', 'GrM');
TagName = evalin('base', 'TagName');

metadata = evalin('base', 'metadata');

ExportData.ds           = metadata{1,1};
ExportData.StartStopInd = metadata{2,1};
ExportData.FitRange     = metadata{3,1};
ExportData.damping      = metadata{4,1};
ExportData.rco          = metadata{7,1};
ExportData.maxCorrRange = metadata{8,1};
ExportData.numMean      = metadata{9,1};
ExportData.CorrSteps    = metadata{10,1};

ExportData.q        = metadata{5,1};
ExportData.r        = metadata{6,1};
ExportData.Sq       = cell2mat(SqM');
ExportData.Gr       = cell2mat(GrM');
ExportData.TagName  = TagName';

%% Export Iq, phiq, Gr Figure
ax4 = handles.axes4;
ax6 = handles.axes6;
ax7 = handles.axes7;

% axes should not be child of panel2, which will be deleted in the figure
set(ax4, 'Parent', gcf);
set(ax6, 'Parent', gcf);
set(ax7, 'Parent', gcf);

% create new figure with isolated axes
A = isolate_axes([ax4, ax6, ax7], 0);

Cs = get(A);
delete(Cs.Children(1));  % delete overlaying panel

set(ax4, 'Parent', handles.Panel2);
set(ax6, 'Parent', handles.Panel2);
set(ax7, 'Parent', handles.Panel2);

% rearranging the axes
set(Cs.Children(7), 'OuterPosition',  [-0.05,   0.6, 1.15, 0.41]);
set(Cs.Children(7), 'XTickLabel',  {''});
set(Cs.Children(7).Title, 'String', ['\textsf{' TagTitle '}']);
set(Cs.Children(5), 'OuterPosition',  [-0.05,   0.3, 1.15, 0.41]);
set(Cs.Children(3), 'OuterPosition',  [-0.05, -0.02, 1.15, 0.41]);

% scaling the window
set(gcf, 'Name', TagTitle);
set(gcf, 'MenuBar', 'figure');
set(gcf, 'OuterPosition', [0.3 0.1 0.45 0.85]);
savefig(gcf, [handles.datpath, '\', char(TagTitle), '.fig'], 'compact');
openfig([handles.datpath, '\', char(TagTitle), '.fig']);
set(gcf, 'Visible', 'on');

% Set Fontsize
set(findall(gcf,'-property','FontSize'), 'FontSize', 14)
% Set linewidth
set(findall(gcf,'-property','Linewidth'), 'Linewidth', 0.8)

delete([handles.datpath, '\', char(TagTitle), '.fig']);

guidata(hObject,handles)

% save figure if checkbox "save" is activated
if handles.save == 1
    extensions = {'fig', 'svg', 'png'};
    rehash
    % create folder name (e.g. _02_03_04_05_06_Export_I_phi_G)
    res_folder_name = '_';
    for fNo = 1:fileNumber
        whitespace = find(isspace(TagName{fNo})==1, 1);
        res_folder_name = [res_folder_name, ...
            char(TagName{fNo}(1:whitespace-1)), ' '];
    end
    res_folder_name = [res_folder_name, 'Export_I_phi_G'];
    res_folder_name = strrep(res_folder_name, ' ', '_');

    % delete latest figure if checkbox "overwrite" is activated
    % and if existing
    rehash
    if handles.overwrite == 1 && ...
            isfolder([handles.datpath, '\', res_folder_name]) == 1

        % find latest figure in folder
        FileInfo = dir([handles.datpath, '\', res_folder_name, '\*.fig']);
        TimeStamp = 0;
        delind = 0;
        for nfile = 1:size(FileInfo)
            FI = FileInfo(nfile,1);
            if FI.datenum > TimeStamp
                TimeStamp = FI.datenum;
                delind = nfile;
            end
        end

        % delete latest figure in folder
        FI = FileInfo(delind,1);
        for k = 1:length(extensions)
            delete([FI.folder, '\', FI.name(1:end-3), extensions{k}]);
        end
        delete([FI.folder, '\ExportData_', FI.name(1:end-3), 'mat']);

        % create dir if not existing
    elseif isfolder([handles.datpath, '\', res_folder_name]) == 0
        mkdir([handles.datpath, '\', res_folder_name])
    end

    % get timestamp
    date = strrep(datestr(now), ':', '-');
    date = strrep(date, ' ', '_');

    % save the figure (02-Apr-2020_16-36-26.fig)
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        date];
    for k = 1:length(extensions)
        saveas(gcf, savename, extensions{k})
    end

    % save ExportData
    savename = [handles.datpath, '\', res_folder_name, '\', ...
        'ExportData_', date,'.mat'];
    save(savename, '-struct', 'ExportData');

    % save Results.xlsx
    copyfile('results_temp.xlsx', [handles.datpath, '\', res_folder_name, '\', ...
        'Results_', date,'.xlsx']);

    disp('Export done.')
end


function plot_hline(val, valE, x1, x2, Color)
% plot reference line
xlims = get(gca, 'xlim');
hold on
plot([x1 x2], [val      val     ],   'k', 'HandleVisibility', 'off', ...
    'Color', Color);
plot([x1 x2], [val+valE val+valE], '--k', 'HandleVisibility', 'off', ...
    'Color', Color);
plot([x1 x2], [val-valE val-valE], '--k', 'HandleVisibility', 'off', ...
    'Color', Color);
hold off


function color = element_color(PSE_Sym)
%% Import VESTA atom colors
filename = 'elements.txt';
formatSpec = '%3f%3s%*6*s%*6*s%*7f%11f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
dataArray{2} = strtrim(dataArray{2});
fclose(fileID);
elements = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','H','VarName6','VarName7','VarName8'});
clearvars filename formatSpec fileID dataArray ans;

Ar  = table2array(elements);
Z   = str2double(Ar(:,1));
El  = Ar(:,2);
RGB = str2double(Ar(:,3:5));

color = [1, 1, 1];
i = 1;
while color == [1, 1, 1]
    if i <= 97 && string(PSE_Sym) == El(i)
        color = [RGB(i,1), RGB(i,2), RGB(i,3)];
    elseif i == 97
        color = [0.1,0.1,0.1];
    end
    i = i + 1;
end