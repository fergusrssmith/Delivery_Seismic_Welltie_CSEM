#!/bin/sh

project_directory=/home/ahglim/disk152/brazil
mincdp=124
inccdp=1
ncdp=1358
minep=2200
incep=1
nep=310
input_filename=hilb_seismic
output_filename=east_bhpintegrate
mincdpevent=1250
maxcdpevent=1481
minepevent=2200
maxepevent=2350
picktop=east_trough.xyz
pickbottom=east_peak.xyz
shifttop=-16
shiftbottom=16
afile=bhpintegrate/east
search=5

sh /home/ahryat/DEdit/east_bhpintegrate_from_editor1a.1.sh
