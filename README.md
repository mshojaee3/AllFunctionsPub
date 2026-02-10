# AllFunctionsPub

**AllFunctionsPub** is a public MATLAB utility package that collects commonly used helper functions for simulation workflows, with a focus on **Abaqus automation**, **data generation**, and **safe directory management**.

The repository is designed to be easily pulled into other projects and used as a MATLAB package (`+AllFunctions`) without manual copying.

---

## 📦 Package Structure

After installation, the package will appear as:




All functions are accessed using the `AllFunctions.` namespace.

---

## 🚀 Installation (Recommended)

Use the following MATLAB script to automatically download the latest versions of the functions into your project:



```matlab

% Root directory = folder of this script (robust)
rootDir = fileparts(mfilename("fullpath"));
addpath(rootDir);   % your helper MATLAB functions are here


% Packages
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
