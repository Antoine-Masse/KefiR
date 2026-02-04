@echo off
setlocal

rem Try to find Rscript in standard locations
set "R_CMD=Rscript"

rem Check if Rscript is in PATH
where Rscript >nul 2>&1
if %ERRORLEVEL% EQU 0 goto Found

rem Check C:\Program Files\R
if exist "C:\Program Files\R" (
    for /f "delims=" %%D in ('dir /b /ad /o-n "C:\Program Files\R"') do (
        if exist "C:\Program Files\R\%%D\bin\Rscript.exe" (
            set "R_CMD=C:\Program Files\R\%%D\bin\Rscript.exe"
            goto Found
        )
    )
)

echo "Rscript not found in PATH or standard locations."
exit /b 1

:Found
echo "Using Rscript at: %R_CMD%"
"%R_CMD%" -e "devtools::check('.', document = FALSE, error_on = 'error')"
