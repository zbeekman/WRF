:: Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

:: Adapted from https://github.com/meta-toolkit/meta/blob/c7019401185cdfa15e1193aad821894c35a83e3f/.appveyor.yml

set THIS_FOLDER=%~dp0

rem Install MSMPI
git clone https://github.com/WRF-CMake/msmpi.git || goto :error
msmpi\MSMpiSetup.exe -unattend || goto :error
msmpi\msmpisdk.msi /passive || goto :error

:: Remember original PATH as refreshenv will replace it with the state from the system registry
set OLDPATH=%PATH%

:: MSMPI addds global vars like MSMPI_INC -- need to propagate these changes to current shell.
:: Use Chocolatey cmd (https://stackoverflow.com/a/44807922).
:: Note: This overrides PATH, hence PATH is updated after this.
call refreshenv.cmd

:: Restore old PATH after refreshenv in case other scripts modified it already
set PATH=%OLDPATH%

:: Setup MinGW shell
set PATH=C:\msys64\usr\bin;%PATH%
set MSYSTEM=MINGW64

rem Install MSYS2 and MinGW64 Packages
bash -lc "pacman --noconfirm -Syu" || goto :error
bash -lc "pacman --noconfirm --needed -S mingw-w64-x86_64-toolchain mingw64/mingw-w64-x86_64-cmake make unzip mingw64/mingw-w64-x86_64-portablexdr" || goto :error
bash -lc "pacman --noconfirm --needed -S mingw-w64-x86_64-libpng mingw-w64-x86_64-libjpeg-turbo mingw-w64-x86_64-jasper" || goto :error
bash -lc "pacman --noconfirm --needed -S mingw-w64-x86_64-hdf5 mingw-w64-x86_64-libtool tar" || goto :error

rem Patch MSMPI
bash -l "%THIS_FOLDER%patch-msmpi.sh" || goto :error
rem Test MSMPI
bash -l "%THIS_FOLDER%test-msmpi.sh" || goto :error
rem Install NetCDF
bash -l "%THIS_FOLDER%install-netcdf.sh" || goto :error

goto :EOF

:error
exit /b %errorlevel%