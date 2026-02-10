function [caseName, params] = getLoadCase(caseNumber)
%GETLOADCASE  Load case from a row in an 8D txt dataset.
% Column order in txtFile:
%   [H11  H22  H12  K11  K22  K12  H31  H32]
%
% Scaling rules (keep sign from txt):
%   H11,H22: multiply by 0.8 if value>=0, else by 0.6
%   H12:     multiply by 0.6
%   K11,K22,K12: multiply by 8 if value>=0, else by 6 
%   H31,H32: multiply by 1.2

% ---- file path (same folder as this function) ----
thisDir = fileparts(mfilename("fullpath"));
txtFile = fullfile(thisDir, "256.txt");
assert(isfile(txtFile), "Cannot find dataset file: %s", txtFile);

% ---- persistent cache ----
persistent Xcache filecache

% Optional: call getLoadCase(0) to reset cache
if caseNumber == 0
    Xcache = [];
    filecache = "";
    caseName = "CACHE_CLEARED";
    params = struct();
    return;
end

% ---- read once ----
if isempty(Xcache) || ~strcmp(filecache, txtFile)
    Xcache = readmatrix(txtFile);
    filecache = txtFile;
    if size(Xcache,2) ~= 8
        error("Expected 8 columns in %s: [H11 H22 H12 K11 K22 K12 H31 H32]. Got %d.", ...
            txtFile, size(Xcache,2));
    end
end

N = size(Xcache,1);
if caseNumber < 1 || caseNumber > N
    error("Invalid caseNumber=%d. Valid range is 1..%d (rows in %s).", caseNumber, N, txtFile);
end

% ---- output defaults ----
params = struct( ...
    'H11_in',1.0,'H22_in',1.0,'H12_in',0.0,'H31_in',0.0,'H32_in',0.0, ...
    'K11_in',0.0,'K22_in',0.0,'K12_in',0.0);

% ---- pick base row ----
x = Xcache(caseNumber, :);

% ---- apply scaling ----
% ---- apply scaling ----
y = zeros(1,8);

% H11, H22: 0.8 if >=0 else 0.6
y(1) = x(1) * ( (x(1) >= 0)*0.8 + (x(1) < 0)*0.6 );  % H11
y(2) = x(2) * ( (x(2) >= 0)*0.8 + (x(2) < 0)*0.6 );  % H22

% H12: always 0.6
y(3) = 0.6 * x(3);

% K11, K22, K12: 8 if >=0 else 6   (keep sign from x)
y(4) = x(4) * ( (x(4) >= 0)*8.0 + (x(4) < 0)*6.0 );  % K11
y(5) = x(5) * ( (x(5) >= 0)*8.0 + (x(5) < 0)*6.0 );  % K22
y(6) = x(6) * ( (x(6) >= 0)*8.0 + (x(6) < 0)*6.0 );  % K12

% H31, H32: always 1.2
y(7) = 1.2 * x(7);
y(8) = 1.2 * x(8);


% ---- map to params ----
addIdentityToH11H22 = true;  % set false if you don't want 1+...
if addIdentityToH11H22
    params.H11_in = 1.0 + y(1);
    params.H22_in = 1.0 + y(2);
else
    params.H11_in = y(1);
    params.H22_in = y(2);
end

params.H12_in = y(3);
params.K11_in = y(4);
params.K22_in = y(5);
params.K12_in = y(6);
params.H31_in = y(7);
params.H32_in = y(8);

% ---- numeric case name/tag ----
caseName = sprintf("RS2_Elem6PerLayer_Case_256_%05d", caseNumber);
end
