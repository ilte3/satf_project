#ifndef __SATF_DEBUG_H__
#define __SATF_DEBUG_H__

class DebugMethod {
  public:
      DebugMethod(const char* class_name, const char* fn_name, int dbg_level, int at_level);
      ~DebugMethod();
      void init();
      
      void log(int at_level, const char* format, ...);
      void set_level(int at_level);
      
  public:
      std::string mClassName;
      std::string mFunctionName;
      int mDbgLevel;
      int mAtLevel;
      
      static int mIndentationLevel;
};

//#define DEBUG

#ifdef DEBUG
  #define _dbg_class_init                   static const char* m_dbg_class; static const int m_dbg_level
  #define _dbg_class_set(class_name, name, dbg_level) const char* class_name::m_dbg_class = name; const int class_name::m_dbg_level = dbg_level
  #define _dbg_method(fn_name, at_level)    DebugMethod dbg_function(m_dbg_class, fn_name, m_dbg_level, at_level)
  #define _dbg_function(fn_name, at_level)  DebugMethod dbg_function("global", fn_name, global_dbg_level, at_level)
  #define _dbg_set_level(at_level)          dbg_function.set_level(at_level)
  #define _dbg(param)                       dbg_function.log param
#else
  #define _dbg_class_init                   
  #define _dbg_class_set(class_name, name, dbg_level) 
  #define _dbg_method(fn_name, at_level)  
  #define _dbg_function(fn_name, at_level)
  #define _dbg_set_level(at_level)
  #define _dbg(param)                              
#endif

#ifdef PROFILING
  #include <unistd.h>
  #define _prof_init    int prof_start
  #define _prof_start   prof_start = getMilliCount()
  #define _prof_add(var)   var += getMilliSpan(prof_start)
#endif

#endif  //__SATF_DEBUG_H__
