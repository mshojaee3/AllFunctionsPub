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
F = ["Mat_1A_ResultantStress_Integral_RVE.m", ...
     "Mat_2A_FEM_2D_nodal_gradient.m", ...
     "Mat_2B_FEM_2D_GaussPoint_gradient.m", ...
     "Mat_3A_FEM_2D_GaussPoint_Integral.m", ...
     "Mat_4A_getLoadCase_PANN_Data_Generation.m", ...
     "Mat_5A_safeCleanRunDir.m", ...
     "Mat_6A_safeCleanRunDirrunAbaqusNoGUI.m", ...
     "Mat_7A_updatePyFromParams.m"];

arrayfun(@(f) websave(fullfile(pkgDir,f), U+f), F, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
