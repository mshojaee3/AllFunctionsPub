function varargout = supabase(action, varargin)
%SUPABASE  All Supabase coordination for the PANN sweep, in one file.
%
% The URL and the key live here and nowhere else. Every machine
% running the sweep shares one case_runs table; before a case is
% run its case_id is looked up, and if a row already exists the
% case is running or finished elsewhere and is skipped.
%
% ACTIONS
% -------
%   [taken, info, ok] = supabase('taken', caseID)
%       taken : true  -> somebody has this case, skip it
%       info  : 'free', 'running on AGAVE', 'completed on EIONE',
%               or the error text when ok is false
%       ok    : false -> the request itself failed. STOP, do not
%               skip: treating an outage as "free" would make
%               every machine run every case.
%
%   ok = supabase('start',  caseID, computerID)
%   ok = supabase('done',   caseID, computerID, message)
%   ok = supabase('fail',   caseID, computerID, message)
%   ok = supabase('delete', caseID)          release a stuck case
%
%   supabase('test')        round-trip check, prints PASS or FAIL
%   supabase('status')      one-line progress summary
%
% EXAMPLE
% -------
%   [taken, info] = supabase('taken', 52);
%   if taken
%       fprintf('SKIP: caseID 52 already %s\n', info);
%   else
%       supabase('start', 52, 'AGAVE');
%       ... run Abaqus ...
%       supabase('done', 52, 'AGAVE', '11/11 frames ok');
%   end
%
% NOTE ON EXACTNESS
% -----------------
% 'taken' followed by 'start' is a check-then-write, not a lock.
% Two machines asking within the same second can both be told
% "free". To make collisions impossible, run once in the Supabase
% SQL editor:
%
%   create unique index case_runs_case_id_uniq
%     on public.case_runs (case_id);
%
% Nothing here needs changing if you add it; the losing insert is
% simply rejected.
%
% Everything below uses curl through system(), not webread/webwrite.

% ============================================================
% CONFIGURATION - the only place these appear
% ============================================================
URL = 'https://mreweybhnsyytsfzxxfj.supabase.co';
KEY = 'sb_publishable_ZtKwPWDEyxQyBy05Kz08oQ_F_agHBZO';

TABLE = 'case_runs2';
% ============================================================

switch lower(char(action))

    case 'taken'
        [t, info, ok] = doTaken(URL, KEY, TABLE, varargin{1});
        varargout = {t, info, ok};

    case 'start'
        ok = doWrite(URL, KEY, TABLE, 'POST', '', struct( ...
                'case_id',     varargin{1}, ...
                'computer_id', char(varargin{2}), ...
                'status',      'running', ...
                'started_at',  utcNow()));
        varargout = {ok};

    case {'done', 'completed'}
        varargout = {doFinish(URL, KEY, TABLE, 'completed', varargin{:})};

    case {'fail', 'failed'}
        varargout = {doFinish(URL, KEY, TABLE, 'failed', varargin{:})};

    case 'delete'
        q  = sprintf('?case_id=eq.%d', varargin{1});
        ok = doWrite(URL, KEY, TABLE, 'DELETE', q, []);
        varargout = {ok};

    case 'test'
        doTest(URL, KEY, TABLE);
        varargout = {};

    case 'status'
        doStatus(URL, KEY, TABLE);
        varargout = {};

    otherwise
        error(['Unknown action "%s". Use: taken, start, done, fail, ' ...
               'delete, test, status.'], char(action));
end
end


% ------------------------------------------------------------
function [taken, info, ok] = doTaken(URL, KEY, TABLE, caseID)
taken = false;  info = '';  ok = false;

url = sprintf('%s/rest/v1/%s?case_id=eq.%d&select=status,computer_id&limit=1', ...
              URL, TABLE, caseID);

cmd = sprintf('curl -sS -H "apikey: %s" -H "Authorization: Bearer %s" "%s"', ...
              KEY, KEY, url);

[st, resp] = system(cmd);
if st ~= 0
    info = sprintf('curl failed (status %d): %s', st, strtrim(resp));
    return;
end

resp = strtrim(resp);
if isempty(resp)
    info = 'empty response from Supabase';
    return;
end

% A result is a JSON array; an error is a JSON object.
if resp(1) ~= '['
    info = sprintf('Supabase error: %s', resp);
    return;
end

ok = true;

if strcmp(resp, '[]')
    info = 'free';
    return;
end

try
    data = jsondecode(resp);
catch ME
    ok = false;
    info = sprintf('could not parse: %s (%s)', resp, ME.message);
    return;
end

if isempty(data)
    info = 'free';
    return;
end

rec = data(1);
stt = 'unknown';   who = 'another computer';
if isfield(rec,'status')      && ~isempty(rec.status);      stt = char(rec.status);      end
if isfield(rec,'computer_id') && ~isempty(rec.computer_id); who = char(rec.computer_id); end

taken = true;
info  = sprintf('%s on %s', stt, who);
end


% ------------------------------------------------------------
function ok = doFinish(URL, KEY, TABLE, newStatus, caseID, computerID, message)
if nargin < 7 || isempty(message); message = ''; end
message = char(message);
if numel(message) > 400; message = [message(1:400) ' ...']; end

q = sprintf('?case_id=eq.%d', caseID);

% result_data is jsonb, so the text is wrapped in an object. If it
% is plain text in your project, pass `message` instead of the struct.
ok = doWrite(URL, KEY, TABLE, 'PATCH', q, struct( ...
        'status',      newStatus, ...
        'computer_id', char(computerID), ...
        'finished_at', utcNow(), ...
        'result_data', struct('message', message)));
end


% ------------------------------------------------------------
function ok = doWrite(URL, KEY, TABLE, verb, query, payloadStruct)
ok = false;

url = sprintf('%s/rest/v1/%s%s', URL, TABLE, query);

dataArg     = '';
payloadFile = '';
if ~isempty(payloadStruct)
    % The JSON goes through a temp file: quoting a brace-and-quote
    % payload inside cmd.exe is fragile, and -d @file behaves the
    % same on Windows, Linux and macOS.
    payloadFile = [tempname '.json'];
    fid = fopen(payloadFile, 'w');
    if fid < 0
        warning('supabase:tmp', 'Could not write payload file.');
        return;
    end
    fwrite(fid, jsonencode(payloadStruct), 'char');
    fclose(fid);
    dataArg = sprintf('-H "Content-Type: application/json" --data-binary @"%s" ', ...
                      payloadFile);
end

bodyFile = [tempname '.txt'];

cmd = sprintf([ ...
    'curl -sS -o "%s" -w "%%{http_code}" -X %s ' ...
    '-H "apikey: %s" -H "Authorization: Bearer %s" ' ...
    '-H "Prefer: return=minimal" %s"%s"'], ...
    bodyFile, verb, KEY, KEY, dataArg, url);

[st, out] = system(cmd);

body = '';
if isfile(bodyFile);              body = strtrim(fileread(bodyFile)); delete(bodyFile); end
if ~isempty(payloadFile) && isfile(payloadFile); delete(payloadFile); end

if st ~= 0
    warning('supabase:curl', 'curl failed (status %d): %s', st, strtrim(out));
    return;
end

code = str2double(regexp(strtrim(out), '\d{3}$', 'match', 'once'));

if any(code == [200 201 204])
    ok = true;
else
    % 409 on a start means another machine inserted first. Harmless.
    warning('supabase:http', '%s returned HTTP %g: %s', verb, code, body);
end
end


% ------------------------------------------------------------
function doTest(URL, KEY, TABLE)
PROBE = 999999;   % far outside the real case range

fprintf('Supabase round-trip test on %s\n\n', URL);

doWrite(URL, KEY, TABLE, 'DELETE', sprintf('?case_id=eq.%d', PROBE), []);

[t1, i1, ok1] = doTaken(URL, KEY, TABLE, PROBE);
fprintf('  1. before start : %s (%s)\n', tf(t1), i1);
if ~ok1
    fprintf(2, '\nFAIL - the check itself failed:\n  %s\n', i1);
    fprintf(2, 'Check the URL, the key, and the RLS policies.\n');
    return;
end

fprintf('  2. writing start row...\n');
doWrite(URL, KEY, TABLE, 'POST', '', struct( ...
    'case_id', PROBE, 'computer_id', 'TEST_PC', ...
    'status', 'running', 'started_at', utcNow()));

[t2, i2] = doTaken(URL, KEY, TABLE, PROBE);
fprintf('  3. after start  : %s (%s)\n', tf(t2), i2);

fprintf('  4. cleaning up...\n');
doWrite(URL, KEY, TABLE, 'DELETE', sprintf('?case_id=eq.%d', PROBE), []);

fprintf('\n');
if ~t1 && t2
    fprintf('PASS - other computers will skip a started case.\n');
else
    fprintf(2, 'FAIL - expected FREE then TAKEN, got %s then %s.\n', tf(t1), tf(t2));
    fprintf(2, 'If step 3 says FREE, the insert is being rejected -\n');
    fprintf(2, 'check the RLS insert policy on %s.\n', TABLE);
end
end


% ------------------------------------------------------------
function doStatus(URL, KEY, TABLE)
url = sprintf('%s/rest/v1/%s?select=case_id,status,computer_id', URL, TABLE);
cmd = sprintf('curl -sS -H "apikey: %s" -H "Authorization: Bearer %s" "%s"', ...
              KEY, KEY, url);
[st, resp] = system(cmd);
resp = strtrim(resp);

if st ~= 0 || isempty(resp) || resp(1) ~= '['
    fprintf(2, 'Could not read %s: %s\n', TABLE, resp);
    return;
end

data = jsondecode(resp);
if isempty(data)
    fprintf('%s is empty - no case has been started yet.\n', TABLE);
    return;
end

if ~iscell(data); data = num2cell(data); end
stt = cellfun(@(r) string(r.status), data);
who = cellfun(@(r) string(r.computer_id), data);

fprintf('%d rows in %s\n', numel(stt), TABLE);
u = unique(stt);
for i = 1:numel(u)
    fprintf('  %-10s %d\n', u(i), sum(stt == u(i)));
end
fprintf('\nby computer:\n');
u = unique(who);
for i = 1:numel(u)
    fprintf('  %-22s %3d done, %d running\n', u(i), ...
            sum(who == u(i) & stt == "completed"), ...
            sum(who == u(i) & stt == "running"));
end
end


% ------------------------------------------------------------
function s = utcNow()
s = char(datetime('now', 'TimeZone', 'UTC', ...
                  'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
end

function s = tf(b)
if b; s = 'TAKEN'; else; s = 'FREE'; end
end
