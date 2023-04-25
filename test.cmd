@echo off
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: This simple batch script helps you to quickly compile your ray-tracing
:: program on Windows.
::
:: Happy hacking! - eric
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: display a message if no filename was given
if [%1] equ [] (
    echo Usage:    .\test {filename}
    echo Example:  .\test scene0.test
    goto :eof
)

:: show the user what the compile command will look like
echo g++ main.cpp std=c++17 -Wall -Iinclude -o main.exe

:: actually compile the program && run it if compilation succeeds
g++ main.cpp -std=c++17 -Wall -Iinclude -o main.exe && main.exe %1 && start test.png
