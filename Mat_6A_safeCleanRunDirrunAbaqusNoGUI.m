function [status, logFile, cmd] = Mat_6A_runAbaqusNoGUI(pyScript, runDir, varargin)
% runAbaqusNoGUI  Run Abaqus/CAE in noGUI mode for an *existing* Python script.
%
% Usage:
%   [status, logFile] = Mat_6A_runAbaqusNoGUI(pyScript, runDir)
%   [status, logFile] = Mat_6A_runAbaqusNoGUI(pyScript, runDir, 'Name',Value,...)
%
% Inputs:
%   pyScript : path to the Python script to run with Abaqus (existing file)
%   runDir   : working directory to run Abaqus in (existing or will be created)
%
% Name-Value options:
%   'AbaqusBat'    : path to abaqus.bat
%   'IntelSetvars' : path to Intel oneAPI setvars.bat
%   'LogFile'      : log file path (default: fullfile(runDir,'log.txt'))
%   'DoIntel'      : true/false, call Intel setvars before Abaqus (default: true)
%   'IntelArch'    : e.g. "intel64" (default: "intel64")
%   'CheckWhere'   : true/false, run "where abaqus" first (default: false)
%   'PrintCmd'     : true/false, display command (default: true)
%
% Outputs:
%   status  : exit code returned by system(cmd) (0 means OK)
%   logFile : path to log file (stdout+stderr captured)
%   cmd     : the command that was executed (useful for debugging)

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
opt.IntelArch    = "intel64";
opt.CheckWhere   = false;
opt.PrintCmd     = true;

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

% ---------- helpers ----------
stripq = @(p) regexprep(char(p), '^"+|"+$', '');  % remove leading/trailing "
q      = @(p) ['"' stripq(p) '"'];                % add exactly one pair of quotes

% ---------- sanitize inputs BEFORE building cmd ----------
runDir          = stripq(runDir);
opt.AbaqusBat    = stripq(opt.AbaqusBat);
opt.IntelSetvars = stripq(opt.IntelSetvars);

% ---------- log file ----------
if (isstring(opt.LogFile) && strlength(opt.LogFile)==0) || (ischar(opt.LogFile) && isempty(opt.LogFile))
    logFile = fullfile(runDir, 'log.txt');
else
    logFile = stripq(opt.LogFile);
end

% ---------- absolute path to python script ----------
pyScriptAbs = pyScript;
try
    pyScriptAbs = char(java.io.File(pyScript).getCanonicalPath());
catch
    % keep as-is if Java not available
end
pyScriptAbs = stripq(pyScriptAbs);

% ---------- optional check ----------
if opt.CheckWhere
    system('where abaqus');
end

% ---------- build command (robust quoting) ----------
if opt.DoIntel
    cmdInner = [ ...
        'cd /d ' q(runDir) ...
        ' && call ' q(opt.IntelSetvars) ' ' char(opt.IntelArch) ...
        ' && ' q(opt.AbaqusBat) ' cae noGUI=' q(pyScriptAbs) ...
        ' > ' q(logFile) ' 2>&1' ...
    ];
else
    cmdInner = [ ...
        'cd /d ' q(runDir) ...
        ' && ' q(opt.AbaqusBat) ' cae noGUI=' q(pyScriptAbs) ...
        ' > ' q(logFile) ' 2>&1' ...
    ];
end

cmd = ['cmd /c "' cmdInner '"'];

if opt.PrintCmd
    disp("CMD = " + string(cmd));
end

% ---------- run ----------
status = system(cmd);

if status ~= 0
    warning('Abaqus returned status=%d. See log: %s', status, logFile);
end

end

