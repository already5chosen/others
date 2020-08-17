cl -O2 -W4 -MD -arch:AVX2 ..\main.cpp ..\chol.cpp ..\chol_innerLoop.c -Fe:s_chol.exe
