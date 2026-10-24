@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0"
if not exist bin mkdir bin
rem Eigen 3.2.91 (C:\eigen) needs _HAS_AUTO_PTR_ETC under C++17 (used by the PDGFa cell)
cl /nologo /O2 /EHsc /openmp /std:c++17 /W1 /utf-8 /wd4244 /I C:\eigen /D_HAS_AUTO_PTR_ETC=1 /DEIGEN_DONT_PARALLELIZE /D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS /D_CRT_SECURE_NO_WARNINGS /Fobin\ /Febin\gisim.exe src\main.cpp src\ICC.cpp src\mCSMC.cpp src\SmoothM.cpp src\pdgfra\CellPart.cpp src\pdgfra\PDGFRaCell.cpp
