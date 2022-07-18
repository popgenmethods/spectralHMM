#!/usr/bin/env python

import os
import subprocess

scriptDir = os.path.dirname(os.path.realpath(__file__))

requiredLibs = ["arpack_combined_all.jar", "JSAP-2.1.jar", "lapack_simple.jar"]
libBase = "spectralHMM_lib"
libDir = os.path.join (scriptDir, libBase)

srcDir = os.path.join (scriptDir, "src")

# check whether the libraries are in place
for lib in requiredLibs:
	if (not os.path.isfile (os.path.join (libDir, lib))):
		print ("[LIBRARY] %s is missing" % lib)
		exit(-1)

# set the classpath right
CLASSPATH = '"%s:%s"' % (os.path.join (libDir, "*"), srcDir)


# path to the main .java file
javaFile = os.path.join (scriptDir, "src/edu/berkeley/spectralHMM/oneD/SelectionHMM.java")

# compile it
javacCmd = "javac -classpath %s %s" % (CLASSPATH, javaFile)
print ("[COMPILING]")
print (javacCmd)
out = subprocess.check_output (javacCmd, shell=True)
if (out not in ["", b'']):
	print ("[ERROR] Error during compiling:")
	print (out)
	exit (-1)

# build the jar
jarFile = os.path.join (scriptDir, "spectralHMM.jar")
jarCmd = "jar cf %s -C %s ." % (jarFile, srcDir)
print ("[CREATING_JAR]")
print (jarCmd)
out = subprocess.check_output (jarCmd, shell=True)
if (out not in ["", b'']):
	print ("[ERROR] Error during creation of jar.")
	print (out)
	exit (-1)

print ("[DONE]")
