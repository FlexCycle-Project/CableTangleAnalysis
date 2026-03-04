# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-build"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/Discregrid"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/tmp"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/extern/Discregrid/src/Ext_Discregrid-stamp${cfgdir}") # cfgdir has leading slash
endif()
