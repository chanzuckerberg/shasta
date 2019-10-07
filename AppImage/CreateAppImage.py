#!/usr/bin/python3

"""
This script creates a Shasta AppImage given a
Shasta dynamic executable. It uses the procedures
described here:
https://docs.appimage.org/packaging-guide/index.html
"""

# Import what we need.
import argparse
import os
import shutil
import stat
import tempfile

# Get the directory where this is script is 
# (normally shasta/AppImage).
# This directory must also contain the skeleton AppDir.
scriptDirectory = os.path.dirname(os.path.abspath(__file__))

# Check that this directory contains the skeleton AppDir.
skeletonAppDir = scriptDirectory + '/AppDir'
if not os.path.isdir(skeletonAppDir):
    raise Exception(
        'The skeleton AppDir was not found. ' +
        'This script should be invoked from its normal location shasta/AppImage.')

# Get the arguments.
parser = argparse.ArgumentParser(description='Create a Shasta AppImage.')
parser.add_argument('shastaInstall', type=str, help='The shasta-install directory.')
shastaInstall = parser.parse_args().shastaInstall

# Work in a temporary directory. 
# The temporary directory gets deleted automatically.
with tempfile.TemporaryDirectory() as tmpDirectory:
    
    # Copy the skeleton AppDir to the temporary directory.
    shutil.copytree(skeletonAppDir, tmpDirectory + '/AppDir')
    
    # Copy the dynamic executable to the AppDir.
    os.mkdir(tmpDirectory + '/AppDir/usr')
    os.mkdir(tmpDirectory + '/AppDir/usr/bin')
    shutil.copy(shastaInstall + '/bin/shastaDynamic', tmpDirectory + '/AppDir/usr/bin/shastaDynamic')
   
    # Copy the shared library to the AppDir.
    os.mkdir(tmpDirectory + '/AppDir/usr/lib')
    shutil.copy(shastaInstall + '/bin/shasta.so', tmpDirectory + '/AppDir/usr/lib/shasta.so')
   
    # Download the linuxdeploy Appimage and make it  executable.
    os.system(
        'wget --quiet ' + 
        'https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage ' + 
        '--directory-prefix ' + tmpDirectory)
    os.chmod(
        tmpDirectory + '/linuxdeploy-x86_64.AppImage',  + 
        stat.S_IRUSR + stat.S_IWUSR + stat.S_IXUSR)

    # Download the appimagetool Appimage and make it  executable.
    os.system(
        'wget --quiet ' + 
        'https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage ' + 
        '--directory-prefix ' + tmpDirectory)
    os.chmod(
        tmpDirectory + '/appimagetool-x86_64.AppImage',  + 
        stat.S_IRUSR + stat.S_IWUSR + stat.S_IXUSR)
    
    # Run linuxdeploy.
    os.system(
        tmpDirectory + 
        '/linuxdeploy-x86_64.AppImage --appdir ' + 
        tmpDirectory + '/AppDir')
        
    # Run appimagetool.
    os.system(
        'ARCH=x86_64 ' + 
        tmpDirectory + 
        '/appimagetool-x86_64.AppImage %s/AppDir %s/shasta-x86_64.AppImage' % (tmpDirectory, tmpDirectory))
        
    # Copy the Shasta AppImage to the desired path.
    shutil.copy(
        '%s/shasta-x86_64.AppImage' % tmpDirectory,
        shastaInstall + '/bin/shasta-x86_64.AppImage')
        
    

