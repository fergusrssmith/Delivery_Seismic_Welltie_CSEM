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


cd ${project_directory}

rm horizons/${afile}.*

bhpread filename=${input_filename} pathlist=int_datasets/${input_filename}.dat keys=ep,cdp keylist=${minepevent}-${maxepevent}:${mincdpevent}-${maxcdpevent} \
| bhploadhdr key=ep,cdp,igc infile=horizons/${picktop} gridtype=R extrap=no bias=${shifttop} \
| bhploadhdr key=ep,cdp,igi infile=horizons/${pickbottom} gridtype=R extrap=no bias=${shiftbottom} \
| sushw key=delrt a=0 b=0 c=0 d=0 j=2147483647 \
| bhpintegrate mode=edge type=trough search=${search} limit=5 \
| bhpstorehdr key=ep,cdp,d2 ofile=horizons/${afile}.area.xyz \
| bhpstorehdr key=ep,cdp,f1 ofile=horizons/${afile}.width.xyz \
| bhpstorehdr key=ep,cdp,f2 ofile=horizons/${afile}.amp.xyz \
| bhpstorehdr key=ep,cdp,igc ofile=horizons/${afile}.p1.xyz \
| bhpstorehdr key=ep,cdp,igi ofile=horizons/${afile}.p2.xyz \
| bhpstorehdr key=ep,cdp,corr ofile=horizons/${afile}.p3.xyz \
| surange

bhpio delete=yes filename=${output_filename} pathlist=int_datasets/${output_filename}.dat

cp int_datasets/default.dat int_datasets/${output_filename}.dat

bhpread filename=${input_filename} pathlist=int_datasets/${input_filename}.dat keys=ep,cdp keylist=*:* \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.area.xyz \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.width.xyz \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.amp.xyz \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.p1.xyz \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.p2.xyz \
| bhphorizon phdr=ep shdr=cdp file=horizons/${afile}.p3.xyz \
| bhpwritecube init=yes key1=cdp,${mincdp},${inccdp},${ncdp} key2=ep,${minep},${incep},${nep} filename=${output_filename} pathlist=int_datasets/${output_filename}.dat horizons=area,width,amp,p1,p2,p3

bhpio delete=yes filename=${output_filename}_transpose_cdp pathlist=int_datasets/${output_filename}_transpose_cdp.dat

bhptranspose filename=${output_filename} pathlist=int_datasets/${output_filename}.dat

cp int_datasets/default.dat int_datasets/${output_filename}_transpose_cdp.dat
