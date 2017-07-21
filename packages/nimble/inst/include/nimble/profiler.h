#ifndef __NIMBLE_PROFILER_H
#define __NIMBLE_PROFILER_H

// This exposes the Google Perftools CPU profiling library.
//   https://github.com/gperftools/gperftools
//   http://gperftools.github.io/gperftools/cpuprofile.html
// These are defined in libprofiler.so, which must be linked in prior to use.
extern "C" int ProfilerStart(const char* fname);
extern "C" void ProfilerStop(void);

#endif //  __NIMBLE_PROFILER_H
