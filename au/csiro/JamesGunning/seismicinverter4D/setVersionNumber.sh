#!/usr/bin/env bash
PKGCLASS=au.csiro.JamesGunning.seismicinverter4D

cat copyright.txt >| VersionNumber.java
echo "package $PKGCLASS;" >> VersionNumber.java
echo "public final class VersionNumber" >> VersionNumber.java
echo "{ public static String versionNumber=\"`date '+%D'`\";}" >> VersionNumber.java

javac -verbose -deprecation VersionNumber.java

