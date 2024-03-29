function [] = visualize(modelfile, datafile)
%VISUALIZE - Animate a 2D mechanism
%
%    VISUALIZE(modelfile, datafile) animates a 2D mechanism, described in
%    the specified ADM file and using data provided in the specified DATA
%    file.
%         modelfile is the name of an ADM file (in JSON format) which is
%                    augmented with visualization information for each
%                    body in the system.
%          datafile  is the name of a ASCII delimited file which contains
%                    one line for each simulation time.  Each line should
%                    contain the value of time, followed by the generalized
%                    coordinates at that time.
%
%     VISUALIZE provides support for turning on/off display of the LRFs and 
%     turning on/off tracking of point trajectories.


% Read model and data

model = loadjson(modelfile);
data = dlmread(datafile);


% Extract the list of bodies from the model
%-------------------------------------------

% Make sure we always work with a cell array.
if iscell(model.bodies)
    blist = model.bodies;
else
    blist = num2cell(model.bodies);
end

nB = length(blist);

% Sanity check: make sure we have data for all bodies in the model
if size(data, 2) ~= 1 + 3*nB
    error('Model/Data mismatch...');
end


% Extract model visualization information
%-----------------------------------------

xl = model.vis.xlim;  % axis xmin/xmax
yl = model.vis.ylim;  % axis ymin/ymax

dx = xl(2)-xl(1);
dy = yl(2)-yl(1);
dd = min(dx, dy);


% Initialize the figure
%-----------------------

hf = figure;
ssize = get(0, 'ScreenSize');
set(hf, 'Position', [60 60 ssize(3)-120 ssize(4)-160]);
hold on
box on
% Draw GRF
quiver([0;0],[0;0],[0.25*dd;0],[0;0.25*dd],'k');
text(0.27*dd, 0, 'x');
text(0, 0.27*dd, 'y');
% Set equal scales in both directions and set axis limits
axis equal
xlim(xl);
ylim(yl);
% Set background color (if specified)
bgColor = [0.9,0.9,0.9];
if isfield(model.vis, 'bgRGB')
    bgColor = model.vis.bgRGB;
end
whitebg(hf, bgColor);
set(hf,'Color',bgColor);
% Hide tick marks and labels
set(gca, 'Xtick', [])
set(gca, 'Ytick', [])

% Create the text object for displaying current time
htexttime = text(xl(1)+dd/20, yl(2)-dd/20, '');
set(htexttime, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Create GUI elements
uicontrol('style', 'PushButton', 'String', 'Quit',...
    'position', [20 20 70 20],...
    'callback', @OnQuit);
uicontrol('style', 'PushButton', 'String', 'Reset',...
    'position', [20 50 70 20],...
    'callback', @OnReset);
uicontrol('style', 'PushButton', 'String', 'Trace On',...
    'position', [20 80 70 20],...
    'callback', @OnTrace);
uicontrol('style', 'PushButton', 'String', 'LRF Off',...
    'position', [20 110 70 20],...
    'callback', @OnLRF);
hslider = uicontrol('style', 'slider', ...
    'Min',0,'Max',1,'Value',0,...
    'Position', [200 20 300 16],...
    'Callback', @OnSlider); 
htextvalue = uicontrol('style', 'text', ...
    'Position', [200, 40 300 20], 'BackgroundColor', bgColor, 'String', '');


% Collect information on all visual objects defined in the model
%----------------------------------------------------------------

% Keep track of the number of each type of visual objects
nS = 0;    % number of shapes
nCG = 0;   % number of CGs
nP = 0;    % number of points

% Loop over all bodies in the model
for iB = 1:nB
    % Current body in list.
    body = blist{iB};
    
    % Skip bodies that do not define visualization information.
    if ~isfield(body, 'vis')
        fprintf('Body %s does not have vis\n', body.name);
        break;
    end
    
    % Append the shapes defined by this body (if any).
    if isfield(body.vis, 'shapes')
        for iS = 1:length(body.vis.shapes)
            % Extract location and size for this shape.
            loc = body.vis.shapes(iS).loc;
            dim = body.vis.shapes(iS).size / 2;
            % Create x and y coordinates for a patch object (expressed in LRF).
            x = [loc(1)-dim(1), loc(1)+dim(1), loc(1)+dim(1), loc(1)-dim(1), loc(1)-dim(1)];
            y = [loc(2)-dim(2), loc(2)-dim(2), loc(2)+dim(2), loc(2)+dim(2), loc(2)-dim(2)];
            % Append a new shape structure.
            nS = nS + 1;
            shapes(nS) = struct('h', [], 'b', body.id, 'color', body.vis.RGB, 'x', x, 'y', y);
        end
    end
    
    % Append this body's CG.
    nCG = nCG + 1;
    cg(nCG) = struct('h', [], 'ha', [], 'b', body.id);
    
    % Append the points defined by this body (if any).
    if isfield(body.vis, 'points')
        for iP = 1:size(body.vis.points, 1)
            x = body.vis.points(iP,1);
            y = body.vis.points(iP,2);
            nP = nP + 1;
            points(nP) = struct('h', [], 'ht', [], 'b', body.id, 'color', body.vis.RGB, 'x', x, 'y', y);
        end
    end
end

% Create all visual objects
%---------------------------
% Note: we do this here, after we collected all objects, to control the
% rendering order: shapes, points, CGs.

for iP = 1:nP
    points(iP).ht = plot(0,0,'k');
    set(points(iP).ht,'Xdata',[],'Ydata',[],'Color',points(iP).color);
end

for iS = 1:nS
    shapes(iS).h = patch(0,0,'b');
    set(shapes(iS).h, 'FaceColor', shapes(iS).color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

for iP = 1:nP
    points(iP).h = plot(0, 0, 'k.');
end

for iCG = 1:nCG
    cg(iCG).h = plot(0, 0, 'ko');
    cg(iCG).ha = quiver([0;0], [0;0], [0.2*dd;0], [0;0.2*dd]);
    set(cg(iCG).ha, 'Color', [0.3 0.3 0.3]);
end


% Process the specified data
t = data(:,1);
nT = length(t);
dt = t(2)-t(1);
q = reshape(data(:,2:end), nT, 3, nB);

% Initial value for the timer period (note: precision = 1 milisecond)
dt = dt - mod(dt, 1e-3);
dt = max(1e-3, dt);

% Set min/max/value for the period slider
set(hslider, 'Min', 1e-3, 'Max', 1, 'Value', dt);
set(htextvalue,'String', sprintf('Refresh: %.3f s', dt));

% Use global variables to make data accessible within callbacks.
global TimerData;
global info;
global iT;
global traj;
global LRF;

info.t = t;
info.q = q;
info.shapes = shapes;
info.points = points;
info.cg = cg;
info.dd = dd;
info.htexttime = htexttime;
info.hslider = hslider;
info.htextvalue = htextvalue;

iT = 0;
traj = false;
LRF = true;

% Create a timer object, then start it
TimerData = timer(...
    'TimerFcn', @OnFrame,...         % callback function (send 'info' as argument)
    'Period',dt,...                  % delay between callback invocations    
    'TasksToExecute',nT,...          % number of calls to callback
    'ExecutionMode','fixedRate',...  % attempt to maintain real-time
    'BusyMode','queue');             % execute callback at next opportunity

start(TimerData);


% ========================================================================
function [] = OnFrame(obj, event)
global info;
global traj;
global iT;

% Unpack info structure.
shapes = info.shapes;
points = info.points;
cg = info.cg;
htexttime = info.htexttime;
t = info.t;
q = info.q;
dd = info.dd;

% Move to next frame.
iT = iT + 1;

if iT > length(q)
    return
end

% Update position of all shapes.
for iS = 1:length(shapes)
    [x,y] = transform2D(shapes(iS).x, shapes(iS).y, q(iT, :, shapes(iS).b));
    set(shapes(iS).h, 'Xdata', x, 'Ydata', y);
end

% Update position of all body-fixed points.
for iP = 1:length(points)
    [x,y] = transform2D(points(iP).x, points(iP).y, q(iT, :, points(iP).b));
    set(points(iP).h, 'Xdata', x, 'Ydata', y);
    if traj
        xx=get(points(iP).ht, 'Xdata');
        yy=get(points(iP).ht, 'Ydata');
        set(points(iP).ht,'Xdata',[xx x], 'Ydata', [yy y]);
    end
end

% Update position of all body CGs and LRFs.
for iCG = 1:length(cg)
    x = q(iT,1,cg(iCG).b);
    y = q(iT,2,cg(iCG).b);
    c = 0.2*dd*cos(q(iT,3,cg(iCG).b));
    s = 0.2*dd*sin(q(iT,3,cg(iCG).b));
    set(cg(iCG).h, 'Xdata', x, 'Ydata', y);
    set(cg(iCG).ha, 'Xdata', [x;x], 'Ydata', [y;y], 'Udata', [c;-s], 'Vdata', [s;c]);
end

% Update the time in the text object
set(htexttime, 'String', sprintf('%4.2f', t(iT)));

% Refresh the figure
drawnow

% =========================================================================
function [xt, yt] = transform2D(x, y, q)
c = cos(q(3));
s = sin(q(3));
xt = q(1) + c * x - s * y;
yt = q(2) + s * x + c * y;

function [xt, yt] = rotate2D(x, y, q)
c = cos(q(3));
s = sin(q(3));
xt = c * x - s * y;
yt = s * x + c * y;


% =========================================================================
% UICONTROLS Callbacks

function [] = OnQuit(obj, event)
% This callback is invoked when the 'Quit' button is pressed.
% Stop and delete the timer object, then close the figure window.
global TimerData;
stop(TimerData);
delete(TimerData);
close(gcf)

function [] = OnTrace(obj, event)
% This callback is invoked to toggle point trajectories.
global traj;
global info;
points = info.points;
for iP = 1:length(points)
    set(points(iP).ht, 'Xdata',[],'Ydata',[]);
end
traj = ~traj;
if traj
    set(obj, 'String', 'Trace Off');
else
    set(obj, 'String', 'Trace On');
end
  
function [] = OnLRF(obj, event)
% This callback is invoked to toggle display of the LRFs.
global LRF;
global info;
cg = info.cg;
LRF = ~LRF;
if LRF
    val = 'on';
    set(obj, 'String', 'LRF Off');
else
    val = 'off';
    set(obj, 'String', 'LRF On');
end
for iCG = 1:length(cg)
   set(cg(iCG).ha, 'Visible', val); 
end

function [] = OnReset(obj, event)
% This callback is invoked when the 'Reset' button is pressed.
% We stop the timer, reset the time counter, and then start the timer.
% Note that we clear all trajectories.
global TimerData;
global info;
global iT;
stop(TimerData);
iT = 0;
points = info.points;
for iP = 1:length(points)
    set(points(iP).ht, 'Xdata',[],'Ydata',[]);
end
start(TimerData);

function [] = OnSlider(obj, event)
% This callback is invoked when the slider value changes.
global TimerData;
global info;

% Read the slider value and adjust to 1 milisecond resolution
dt = get(info.hslider,'Value');
dt = dt - mod(dt,.001);

% Set slider readout to show its value
set(info.htextvalue,'String', sprintf('Refresh: %.3f s', dt));

% If timer is on, stop it, reset the period, and start it again.
if strcmp(get(TimerData, 'Running'), 'on')
    stop(TimerData);
    set(TimerData,'Period',dt)
    start(TimerData)
else               % If timer is stopped, reset its period only.
    set(TimerData,'Period',dt)
end
