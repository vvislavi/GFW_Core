#!/bin/bash
function AddToCMakeFile() { #1 - Line to add, 2 - before which one, 3 - file name
    #Check if the line to add already exists
    lex=`grep ${1} ${3}`
    test ! -z "${lex}" && return
    #Find where we have to return it
    lex=`grep ${2} ${3}`
    test -z "${lex}" && return
    sed -i '' -E "s/([[:space:]]+)${2}/\1${1}NLPlaceHolder\1${2}/g" ${3}
    sed -i '' -E 's/NLPlaceHolder/\'$'\n/g' ${3}
    unset lex
}
function FixTask() {
  #Remove power arrays
  sed -i '' "/int pows/d" $1
  #Remove power array arguments from the file, too
  sed -i '' "s/ 7, powsFull,//g" $1
  sed -i '' "s/ 7, pows,//g" $1
  #Change from Head.Data() to Head.c_str
  sed -i '' "s/Head.Data()/Head.c_str()/g" $1
  #Change from Re() and Im() to real() and im $1ag()
  sed -i '' "s/Re()/real()/g" $1
  sed -i '' "s/Im()/imag()/g" $1
  #Change TComplex to complex<double>
  sed -i '' "s/TComplex/complex\<double\>/g" $1
  #Add CreateRegions() after the last corrconfig has been added
  lastCFG='corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP {2 3} refN {-2 -3}", "ChSC234", kFALSE));'
  lastLine=`sed -n "/${lastCFG}/=" $1`
  if [ ! -z ${lastLine} ]
  then
    ((lastLine++))
    sed -i '' "${lastLine} i\\
    fGFW->CreateRegions();\\
    " $1
  fi
  #Also, dqTask has explicitly preset DisableOverlap flag. Remove that
  sed -i '' "/bool DisableOverlap/d" $1
  sed -i '' "s/, DisableOverlap//g" $1
  unset lastCFG
  unset lastLine
}
#First, check if the relevant directory is set and if we can figure out the path
if [ -z $1 ]
then
  echo "Please specify the O2Physics directory!"
  return
fi
tarDir=$1
if [ ! -d ${1}/PWGCF ]
then
  if [ -d ${1}/O2Physics/PWGCF ]
  then
    tarDir=${1}/O2Physics
  else
    echo "Could not find O2Physics in the directory specified!"
    return
  fi
fi
#copy the relevant files
cp GFW*.cxx ${tarDir}/PWGCF/GenericFramework/
cp GFW*.h ${tarDir}/PWGCF/GenericFramework/
AddToCMakeFile GFWPowerArray.cxx GFW.cxx ${tarDir}/PWGCF/GenericFramework/CMakeLists.txt
AddToCMakeFile GFWPowerArray.h GFW.h ${tarDir}/PWGCF/GenericFramework/CMakeLists.txt
FixTask ${tarDir}/PWGCF/Tasks/flowGenericFramework.cxx
FixTask ${tarDir}/PWGDQ/Core/VarManager.h
FixTask ${tarDir}/PWGDQ/Tasks/dqFlow.cxx
echo "All done! To stage all changes for commit, please run:"
echo "cd ${tarDir} && git add PWGCF/GenericFramework/GFWPowerArray.cxx PWGCF/GenericFramework/GFWPowerArray.h PWGCF/GenericFramework/GFW.cxx PWGCF/GenericFramework/GFW.h PWGCF/GenericFramework/GFWCumulant.cxx PWGCF/GenericFramework/GFWCumulant.h PWGCF/GenericFramework/CMakeLists.txt PWGCF/Tasks/flowGenericFramework.cxx PWGDQ/Core/VarManager.h PWGDQ/Tasks/dqFlow.cxx"
