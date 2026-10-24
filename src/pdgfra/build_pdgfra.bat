@echo off
rem Build pdgfra_test.exe (headless PDGFRa cell test) with MSVC 2022.
rem C:\eigen is Eigen 3.2.91 (pre-3.3), which uses std::binder1st/binder2nd removed
rem in C++17; _HAS_AUTO_PTR_ETC=1 restores them in the MSVC STL.
rem Invoke from PowerShell:  & cmd.exe /c "`"<fullpath>\build_pdgfra.bat`""
setlocal
cd /d "%~dp0"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
if not exist "..\..\bin" mkdir "..\..\bin"
if not exist "obj" mkdir "obj"
cl /nologo /O2 /EHsc /openmp /std:c++17 /utf-8 /W3 /wd4244 /I C:\eigen /D_HAS_AUTO_PTR_ETC=1 /D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS /D_CRT_SECURE_NO_WARNINGS /Fo"obj\\" ^
   CellPart.cpp PDGFRaCell.cpp pdgfra_test.cpp /Fe"..\..\bin\pdgfra_test.exe"
set ERR=%ERRORLEVEL%
if %ERR% neq 0 echo BUILD FAILED (%ERR%)
exit /b %ERR%
