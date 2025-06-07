@echo off
setlocal EnableDelayedExpansion

REM This script builds a Python package using setuptools and creates wheel distributions

REM Clean up previous builds
echo Cleaning up previous builds...
if exist "build" rd /s /q "build"
if exist "dist" rd /s /q "dist"
if exist "microaggregation.egg-info" rd /s /q "microaggregation.egg-info"

REM Loop through Python environments
for %%v in (py37 py38 py39 py310 py311 py312) do (
    echo.
    echo ========================================
    echo Building for Python environment: %%v
    echo ========================================
    
    call conda activate %%v
    if !errorlevel! neq 0 (
        echo Error: Failed to activate %%v environment
        goto :error_exit
    )

    REM Ensure build tools are present and reasonably up-to-date in each environment
    echo Installing/upgrading build tools for %%v...
    python -m pip install --upgrade pip build setuptools wheel pybind11 numpy
    if !errorlevel! neq 0 (
        echo Error: Failed to install/upgrade build tools for %%v
        call conda deactivate
        goto :error_exit
    )

    REM Use python -m build for sdist and wheel
    echo Building sdist and wheel for %%v...
    python -m build
    if !errorlevel! neq 0 (
        echo Error: Build failed for %%v using 'python -m build'
        call conda deactivate
        goto :error_exit
    )
    
    call conda deactivate
)

REM Upload the distributions to PyPI
echo.
echo ========================================
echo Uploading distributions to PyPI...
echo ========================================

REM Use a consistent, preferably newer, environment for uploading
call conda activate py311 
if !errorlevel! neq 0 (
    echo Error: Failed to activate py311 for uploading
    goto :error_exit
)

echo Installing/upgrading twine...
python -m pip install --upgrade twine
if !errorlevel! neq 0 (
    echo Error: Failed to install/upgrade twine
    call conda deactivate
    goto :error_exit
)

echo Running twine upload...
twine upload dist/* --verbose
if !errorlevel! neq 0 (
    echo Error: Upload to PyPI failed
    call conda deactivate
    goto :error_exit
)

call conda deactivate

echo.
echo ========================================
echo All builds and upload completed successfully!
echo ========================================
goto :eof

:error_exit
echo Build or upload process failed.
exit /b 1

:eof
