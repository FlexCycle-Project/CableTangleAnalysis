# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitclone-lastrun.txt" AND EXISTS "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitinfo.txt" AND
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitclone-lastrun.txt" IS_NEWER_THAN "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://github.com/InteractiveComputerGraphics/Discregrid.git" "Ext_Discregrid"
    WORKING_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/InteractiveComputerGraphics/Discregrid.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"
          checkout "ddf20dc0480874bf02e0bdc6ded76c1f101b17fb" --
  WORKING_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'ddf20dc0480874bf02e0bdc6ded76c1f101b17fb'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitinfo.txt" "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/Ext_Discregrid-gitclone-lastrun.txt'")
endif()
