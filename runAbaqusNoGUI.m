% runDir  = "C:\myproj\AbaqusSimulation";
% pyOut   = fullfile(runDir, "main_updated.py");   % already created by your update function

% [status, logFile] = runAbaqusNoGUI(pyOut, runDir, ...
    % 'AbaqusBat', "C:\SIMULIA\Commands\abaqus.bat", ...
    % 'IntelSetvars', "C:\Program Files (x86)\Intel\oneAPI\setvars.bat", ...
    % 'DoIntel', true);


function [status, logFile] = runAbaqusNoGUI(pyScript, runDir, varargin)
% runAbaqusNoGUI  Run Abaqus/CAE in noGUI mode for an *existing* Python script.
%
% Usage:
%   [status, logFile] = runAbaqusNoGUI(pyScript, runDir)
%   [status, logFile] = runAbaqusNoGUI(pyScript, runDir, 'Name',Value,...)
%
% Inputs:
%   pyScript : path to the Python script to run with Abaqus (existing file)
%   runDir   : working directory to run Abaqus in (existing or will be created)
%
% Name-Value options:
%   'AbaqusBat'    : path to abaqus.bat (default: "C:\SIMULIA\Commands\abaqus.bat")
%   'IntelSetvars' : path to Intel oneAPI setvars.bat (default: "C:\Program Files (x86)\Intel\oneAPI\setvars.bat")
%   'LogFile'      : log file path (default: fullfile(runDir,'log.txt'))
%   'DoIntel'      : true/false, call Intel setvars before Abaqus (default: true)
%   'CheckWhere'   : true/false, run "where abaqus" first (default: false)
%
% Outputs:
%   status  : exit code returned by system(cmd) (0 means OK)
%   logFile : path to log file (stdout+stderr captured)

% ---------- validate ----------
assert(ischar(pyScript) || isstring(pyScript), 'pyScript must be char/string.');
assert(ischar(runDir)   || isstring(runDir),   'runDir must be char/string.');

pyScript = char(pyScript);
runDir   = char(runDir);

assert(isfile(pyScript), 'Python script not found: %s', pyScript);

if ~exist(runDir, 'dir')
    mkdir(runDir);
end

% ---------- defaults ----------
opt.AbaqusBat    = "C:\SIMULIA\Commands\abaqus.bat";
opt.IntelSetvars = "C:\Program Files (x86)\Intel\oneAPI\setvars.bat";
opt.LogFile      = "";
opt.DoIntel      = true;
opt.CheckWhere   = false;

% ---------- parse name-value ----------
if mod(numel(varargin),2) ~= 0
    error('Optional inputs must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    val  = varargin{i+1};
    if ~isfield(opt, name)
        error('Unknown option "%s".', name);
    end
    opt.(name) = val;
end

if strlength(opt.LogFile) == 0
    logFile = fullfile(runDir, 'log.txt');
else
    logFile = char(opt.LogFile);
end

% Use absolute path for the script in the command (safer)
pyScriptAbs = pyScript;
try
    pyScriptAbs = char(java.io.File(pyScript).getCanonicalPath());
catch
    % if java not available for some reason, keep original
end

status = 0;

if opt.CheckWhere
    system('where abaqus');
end

% ---------- build command ----------
if opt.DoIntel
    % IMPORTANT: escape quotes for Windows cmd properly using "" inside "...".
    cmd = sprintf( ...
        'cmd /c "cd /d ""%s"" && call ""%s"" intel64 && ""%s"" cae noGUI=""%s"" > ""%s"" 2>&1"', ...
        runDir, opt.IntelSetvars, opt.AbaqusBat, pyScriptAbs, logFile);
else
    cmd = sprintf( ...
        'cmd /c "cd /d ""%s"" && ""%s"" cae noGUI=""%s"" > ""%s"" 2>&1"', ...
        runDir, opt.AbaqusBat, pyScriptAbs, logFile);
end

% ---------- run ----------
status = system(cmd);

if status ~= 0
    warning('Abaqus returned status=%d. See log: %s', status, logFile);
end

end
