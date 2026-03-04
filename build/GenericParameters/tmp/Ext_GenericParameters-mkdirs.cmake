# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters-build"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/tmp"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters-stamp"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src"
  "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/afs/informatik.uni-goettingen.de/user/s/sezer.kemalalp/Desktop/PositionBasedDynamics-master/PositionBasedDynamics-master/build/GenericParameters/src/Ext_GenericParameters-stamp${cfgdir}") # cfgdir has leading slash
endif()
