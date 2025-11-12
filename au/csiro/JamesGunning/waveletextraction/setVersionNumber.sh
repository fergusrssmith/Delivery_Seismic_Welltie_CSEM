#!/usr/bin/env bash
PKGCLASS=au.csiro.JamesGunning.waveletextraction

cat copyright.txt >| VersionNumber.java
echo "package $PKGCLASS;" >> VersionNumber.java
echo "public final class VersionNumber" >> VersionNumber.java
echo "{ public static String versionNumber=\"`date '+%D'`\";}" >> VersionNumber.java

javac VersionNumber.java

