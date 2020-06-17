Windows variant of Travis Downs's Hardware Store Elimination test.
See https://travisdowns.github.io/blog/2020/05/13/intel-zero-opt.html#our-first-benchmark for details.

The test is built in two variants: travis_test and travis_test_wlp.
The first variant uses regular memory allocation which on current Windows ends up in 4KB pages.
The second variant uses VirtualAlloc() with MEM_LARGE_PAGES flag.

The benchmark has to be compiled with 64-bit clang9 or 10. Other compilers are not tested, so results are not guaranteed to be compatible.
Or should I say, almost guaranteed to be incompatible?

travis_test_wlp is not easy to run, it has to be run from account that has SeLockMemoryPrivilege enabled.
Enabling is done by Local Group Policy Editor (gpedit) that you have to run as administrator.
Location: Computer Configuration\Windows Setting\Security Settings\Local Policies\User Rights Assignment
Key Name: Lock Pages in memory
Click on the key and in open dialog add the account from which you want to run a test.
The settings come into effect after next logon.
I had better success with *regular user account*, *not* admin.

Results:
cfl1-lp.txt - Xeon E-2176G @4.25 GHz, 2xDDR4-2666, large pages
cfl1-sp.txt - Xeon E-2176G @4.25 GHz, 2xDDR4-2666, default pages
