@echo off
setlocal EnableDelayedExpansion

REM This script builds a Python package using setuptools and creates wheel distributions

REM Clean up previous builds
echo Cleaning up previous builds...
if exist "build" rd /s /q "build"
if exist "dist" rd /s /q "dist"
if exist "*.egg-info" rd /s /q "*.egg-info"

REM Loop through Python environments
for %%v in (py37 py38 py39 py310 py311 py312) do (
    echo.
    echo ========================================
    echo Building for Python environment: %%v
    echo ========================================
    
    call conda activate %%v
    if !errorlevel! neq 0 (
        echo Error: Failed to activate %%v environment
        exit /b !errorlevel!
    )

    python setup.py build
    if !errorlevel! neq 0 (
        echo Error: Build failed for %%v
        call conda deactivate
        exit /b !errorlevel!
    )

    python setup.py bdist_wheel
    if !errorlevel! neq 0 (
        echo Error: Wheel creation failed for %%v
        call conda deactivate
        exit /b !errorlevel!
    )
    
    call conda deactivate
)

REM Upload the wheel to PyPI
echo.
echo ========================================
echo Uploading wheels to PyPI...
echo ========================================

call conda activate py311
twine upload dist/* --verbose
if !errorlevel! neq 0 (
    echo Error: Upload to PyPI failed
    call conda deactivate
    exit /b !errorlevel!
)

call conda deactivate

echo.
echo ========================================
echo All builds completed successfully!
echo ========================================
exit /b 0
