@echo off

:: copy repo ssh keys to working dir
IF NOT EXIST "%USERPROFILE%\.ssh" (
	MD "%USERPROFILE%\.ssh"
)
IF EXIST ".ssh\id_rsa" (
	COPY ".ssh\id_rsa" "%USERPROFILE%\.ssh\id_rsa"
)
IF EXIST ".ssh\id_rsa.pub" (
	COPY ".ssh\id_rsa.pub" "%USERPROFILE%\.ssh\id_rsa.pub"
)

:: start matlab project
matlab.exe -r "setenv('HOME', [getenv('USERPROFILE')]); open('GDPPEMELSystem.prj')"