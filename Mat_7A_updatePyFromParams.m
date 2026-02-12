% clc
% clear
% params.H11 = 0.1;
% params.H12 = 0.0;
% params.H22 = 0.0;
% params.H21 = 0.0;
% rootDir = fileparts(mfilename("fullpath"));
% runDir = fullfile(rootDir, "AbaqusSimulation");
% if ~exist(runDir,'dir')
%     mkdir(runDir);
% end
% srcPy = fullfile(rootDir, "main.py");
% AllFunctions.Mat_7A_updatePyFromParams(params, srcPy, fullfile(runDir,"main_updated.py"))
function outPy = Mat_7A_updatePyFromParams(params, pyFile, outPy)

% ----------------- validate -----------------
assert(isstruct(params), 'params must be a struct.');
assert(ischar(pyFile) || isstring(pyFile), 'pyFile must be char/string.');

pyFile = char(pyFile);
assert(isfile(pyFile), 'Cannot find python file: %s', pyFile);

[pyDir, baseName, ext] = fileparts(pyFile);
assert(strcmpi(ext, '.py'), 'pyFile must be a .py file.');

if nargin < 3 || isempty(outPy)
    outPy = fullfile(pyDir, [baseName '_updated.py']);
else
    outPy = char(outPy);
    outDir = fileparts(outPy);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% ----------------- read + split into lines -----------------
txt   = fileread(pyFile);
lines = regexp(txt, '\r\n|\n|\r', 'split');

names = fieldnames(params);

for k = 1:numel(names)
    varName = names{k};
    pyLit   = matlab_to_python_literal(params.(varName));

    safeVar = regexptranslate('escape', varName);

    % match: optional leading whitespace, varName, optional spaces, '=', anything
    patLine = ['^\s*' safeVar '\s*=.*$'];

    found = false;
    for i = 1:numel(lines)
        if ~isempty(regexp(lines{i}, patLine, 'once'))
            % Replace RHS but keep "varName = " prefix
            lines{i} = regexprep(lines{i}, ...
                ['(^\s*' safeVar '\s*=\s*).*$'], ...
                ['$1' pyLit]);
            found = true;
            break;
        end
    end

    if ~found
        if ~isempty(lines) && ~isempty(lines{end})
            lines{end+1} = '';
        end
        lines{end+1} = sprintf('%s = %s', varName, pyLit);
    end
end

outTxt = strjoin(lines, newline);

% ----------------- enforce UTF-8 coding header -----------------
if isempty(regexp(outTxt,'coding[:=]\s*[-\w.]+','once'))
    outTxt = sprintf('# -*- coding: utf-8 -*-\n%s', outTxt);
end

% ----------------- write output -----------------
fid = fopen(outPy, 'w');
assert(fid ~= -1, 'Could not open output file for writing: %s', outPy);
fwrite(fid, outTxt, 'char');
fclose(fid);

end % ===== end main function =====


% ======================================================================
% Local helpers (still ONE file)
% ======================================================================

function lit = matlab_to_python_literal(v)
if islogical(v) && isscalar(v)
    lit = tern(v, 'True', 'False');

elseif isnumeric(v) && isscalar(v)
    lit = sprintf('%.15g', v);

elseif ischar(v) || isstring(v)
    s = char(v);
    s = strrep(s, '\', '\\');
    s = strrep(s, '''', '\''');           % escape single quote for Python string
    lit = ['''' s ''''];

elseif isnumeric(v) && isvector(v)
    elems = arrayfun(@(x) sprintf('%.15g', x), v(:).', 'UniformOutput', false);
    lit = ['[' strjoin(elems, ', ') ']'];

elseif iscellstr(v) || (iscell(v) && all(cellfun(@(x) ischar(x) || isstring(x), v)))
    elems = cellfun(@(x) ['''' strrep(strrep(char(x), '\','\\'),'''','\''') ''''], ...
        v, 'UniformOutput', false);
    lit = ['[' strjoin(elems, ', ') ']'];

else
    error('Unsupported value type for automatic Python literal conversion.');
end
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

