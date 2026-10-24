// Stub replacing the MFC precompiled header so the original cell models build headless.
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
typedef int BOOL;
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
