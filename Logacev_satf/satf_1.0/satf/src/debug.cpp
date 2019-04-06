#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include "debug.h"


int DebugMethod::mIndentationLevel = 0;

DebugMethod::DebugMethod(const char* class_name, const char* fn_name, int dbg_level, int at_level) {
    mDbgLevel = dbg_level;
    mAtLevel = at_level;
    mClassName = class_name;
    mFunctionName = fn_name;
    log(mAtLevel, fn_name);
    mIndentationLevel++;
}

DebugMethod::~DebugMethod() {
    mIndentationLevel--;
    log(mAtLevel, "--");
}

void DebugMethod::set_level(int at_level){
      mAtLevel = at_level;
}

void DebugMethod::log(int at_level, const char * format, ...) 
{
    if(at_level+mAtLevel > mDbgLevel)
      return;
      
    for(int i=0; i < mIndentationLevel; i++)
      putc('\t', stdout);
    
    char buffer[1000];
    va_list args;
    va_start (args, format);
    vsprintf (buffer, format, args);
    printf(" %s", buffer);
    va_end (args);
    printf(" (%s)\n", mClassName.c_str());
}
