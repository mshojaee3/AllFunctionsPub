# AllFunctionsPub

%Packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkgDir  = fullfile(rootDir, "+AllFunctions");
if ~exist(pkgDir,"dir"); mkdir(pkgDir); end

U = "https://raw.githubusercontent.com/mshojaee3/AllFunctionsPub/main/";
F = ["updatePyFromParams.m","runAbaqusNoGUI.m", ...
     "getLoadCase_PANN_Data_Generation.m", ...
     "safeCleanRunDir.m", ...
     "ShellResultantStress_Based_Volume_Integratin_of_RVE_2_10_2026.m"];
arrayfun(@(f) websave(fullfile(pkgDir,f), U+f), F, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
