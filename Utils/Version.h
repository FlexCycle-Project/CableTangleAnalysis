#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "87af141cc03487e353f748af66dc3c717eeb0915"
#define GIT_REFSPEC "refs/heads/linux"
#define GIT_LOCAL_STATUS "CLEAN"

#define PBD_VERSION "2.2.2"

#ifdef DL_OUTPUT

#endif

#endif
