#!/bin/bash

Help()
{
   echo "wrapper for optiso"
   echo
   echo "Syntax: optrun.sh [-d|-p|-c]"
   echo "d  directory with config.json files for individual genes"
   echo "p  path to apc program to use"
   echo "c  number of cpus to use"
   echo
}

while getopts ":hv" option; do
    echo $option 
    case $option in
        h)
            Help
            exit;;
        \?)
            echo "Error: Invalid option"
            exit;;
    esac
done

# defaults
Name = "world"

echo "hell $Name"