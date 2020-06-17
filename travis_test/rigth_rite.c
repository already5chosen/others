#define UNICODE
#define _UNICODE
#include <windows.h>
#include <tchar.h>
#include <stdio.h>

// #define BUF_SIZE 65536

// TCHAR szName[]=TEXT("LARGEPAGE");
// typedef int (*GETLARGEPAGEMINIMUM)(void);

void DisplayError(const wchar_t* pszAPI, DWORD dwError)
{
  LPVOID lpvMessageBuffer;
  FormatMessage(
    FORMAT_MESSAGE_ALLOCATE_BUFFER |
    FORMAT_MESSAGE_FROM_SYSTEM |
    FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL, dwError,
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
    (LPTSTR)&lpvMessageBuffer, 0, NULL);

  //... now display this string
  _tprintf(TEXT("%s() failed with code %d. %s"), pszAPI, dwError, lpvMessageBuffer);

  // Free the buffer allocated by the system
  LocalFree(lpvMessageBuffer);
}

int Privilege(const wchar_t* pszPrivilege, BOOL bEnable)
{
  // open process token
  HANDLE hToken;
  if (!OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &hToken)) {
    DisplayError(TEXT("OpenProcessToken"), GetLastError());
    return 0;
  }

  // get the luid
  TOKEN_PRIVILEGES tp;
  if (!LookupPrivilegeValue(NULL, pszPrivilege, &tp.Privileges[0].Luid)) {
    DisplayError(TEXT("LookupPrivilegeValue"), GetLastError());
    return 0;
  }

  // enable or disable privilege
  tp.PrivilegeCount = 1;
  tp.Privileges[0].Attributes = bEnable ? SE_PRIVILEGE_ENABLED : 0;
  BOOL status = AdjustTokenPrivileges(hToken, FALSE, &tp, 0, (PTOKEN_PRIVILEGES)NULL, 0);
  // It is possible for AdjustTokenPrivileges to return TRUE and still not succeed.
  // So always check for the last error value.
  DWORD error = GetLastError();
  if (!status || (error != ERROR_SUCCESS)) {
    DisplayError(TEXT("AdjustTokenPrivileges"), GetLastError());
    return 0;
  }

  // close the handle
  if (!CloseHandle(hToken)) {
    DisplayError(TEXT("CloseHandle"), GetLastError());
    return 0;
  }

  return 1; // success
}

int RightRite(void)
{
  return Privilege(TEXT("SeLockMemoryPrivilege"), TRUE);
}
