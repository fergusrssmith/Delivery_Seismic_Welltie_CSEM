#!/usr/bin/env bash
PKGCLASS=au.csiro.JamesGunning.CSEM.OneDinversion

cat copyright.txt >| VersionNumber.java
echo "package $PKGCLASS;" >> VersionNumber.java
echo "public final class VersionNumber" >> VersionNumber.java
echo "{ public static String versionNumber=\"`date '+%D'`\";}" >> VersionNumber.java

javac -g -O -verbose -deprecation VersionNumber.java

