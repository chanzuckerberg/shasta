name: Build

on: [push]

jobs:
  build:

    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@master
    - name: InstallPrerequisites
      run: ./scripts/InstallPrerequisites.sh
    - name: Build
      run: (mkdir shasta-build; cd shasta-build; cmake ..; make -j all)
   
