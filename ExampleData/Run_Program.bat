@echo off

set ExePath=PeptideProphetRunner.exe

if exist %ExePath% goto DoWork
if exist ..\%ExePath% set ExePath=..\%ExePath% && goto DoWork
if exist ..\bin\%ExePath% set ExePath=..\bin\%ExePath% && goto DoWork
if exist ..\PeptideProphetRunner\bin\%ExePath% set ExePath=..\PeptideProphetRunner\bin\%ExePath% && goto DoWork
if exist ..\PeptideProphetRunner\bin\Debug\%ExePath% set ExePath=..\PeptideProphetRunner\bin\debug\%ExePath% && goto DoWork

echo Executable not found: %ExePath%
goto Done

:DoWork
echo.
echo Procesing with %ExePath%
echo.

%ExePath% QC_Shew_12_02-3_250ng_1pt9um_75um_Frodo_16Feb13_syn.txt

:Done

pause
