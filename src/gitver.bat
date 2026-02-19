@echo off
setlocal

REM %~dp0 is this script's directory (Windows path, ends with backslash)
REM ProjectDir points to the .vcxproj directory; use that as the base.
set "PROJDIR=%~1"
if "%PROJDIR%"=="" set "PROJDIR=%CD%"

REM Convert Windows path to WSL path and run the script
wsl bash -lc "cd \"$(wslpath -a '%PROJDIR%')\" && ./gitver.bash"
exit /b %ERRORLEVEL%