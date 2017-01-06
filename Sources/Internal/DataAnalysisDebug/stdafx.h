// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <Windows.h>
#include <ppl.h>

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#include <string>
#include <stdio.h>

#include "../../External/ExternalLib.h"
#include "../Utils/UtilsLib.h"
#include "../DataAnalysis/CommonTypes.h"

using namespace DataAnalysis;
using namespace DataAnalysis::Utils;

#include "../Transformations/TransformationsLib.h"
using namespace DataAnalysis::Transformations;

// TODO: reference additional headers your program requires here
