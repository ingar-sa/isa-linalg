@echo off

set BuildFolder="%CD%"\build_windows\

REM set CommonCPPCompilerFlags=-MT -nologo -Gm- -GR- -EHa -Od -Oi -W4 -wd4201 -wd4100 -wd4189 -wd4505 -Z7 /std:c++20 /Fe%BuildFolder% /Fo%BuildFolder%

set CommonCPPCompilerFlags=-MT -nologo -Gm- -GR- -EHa -O2 -Oi -W4 -wd4201 -wd4100 -wd4189 -wd4505 -Z7 /std:c++20 /Fe%BuildFolder% /Fo%BuildFolder%

REM setCommonCCompilerFlags=-MT -nologo -Gm- -GR- -EHa -Od -Oi -W4 -wd4201 -wd4100 -wd4189 -wd4505 -Z7 /Fe%BuildFolder% /Fo%BuildFolder%
set CommonCCompilerFlags= -nologo -MD -Os /Fe%BuildFolder% /Fo%BuildFolder%

For /R %%G in (*.cpp) do cl %CommonCPPCompilerFlags% "%%G"
For /R %%G in (*.c) do cl %CommonCCompilerFlags% "%%G"
