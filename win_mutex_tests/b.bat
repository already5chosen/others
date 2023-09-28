cl -W4 -O2 -MD -nologo mtx_test.c -Fe: crit_sec_test
cl -W4 -O2 -MD -nologo mtx_test.c -Fe: srw_lock_test -DUSE_SRWL
