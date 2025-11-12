#!/usr/bin/env bash
#
# This script converts old XML formats to new ones
#
# Accrued changes are
#
# (1) signal_to_noise -> noise_rms
# (2) layer sequence is now -> top_layer, middle_layer, middle_layer....,bottom_layer

# Check usage

if [ $# -eq 2 ];
then
    if [ $1 == "-u" ];
	then
		# undo conversion
		undo=1
		file=$2
	else
		echo "Usage: $0 [-u] xmlfile"
		echo "supply -u to undo a conversion"
		exit 0
	fi
elif [ $# -eq 1 ];
then
    undo=0
	file=$1
else
	echo "Usage: $0 [-u] xmlfile"
	echo "supply -u to undo a conversion"
	exit 0
fi

# now do the conversion

if [ $undo -eq 1 ];
then
	# Undo is simple  
	# noise_rms -> signal_to_noise
	#
	echo "translating $file to old format"
	(awk '{gsub(/noise_rms/,"signal_to_noise")} \
           /<\/?top_layer/{gsub(/top_layer/,"layer")} \
           /<\/?middle_layer/{gsub(/middle_layer/,"layer")} \
		  {print $0}' \
	   $file > $file.tmp \
	&& mv -f $file.tmp $file ) || echo "Failed translation in: $file"
	exit 0
else
	# forward translation
	#
	# Do signal_to_noise -> noise_rms
	echo "translating $file to new format"
	(awk '/signal_to_noise/{gsub(/signal_to_noise/,"noise_rms")} \
          /<layer/{l++; if (l==1) {gsub(/layer/,"top_layer")} \
                        else {gsub(/layer/,"middle_layer")}} \
          /<\/layer/{if (l==1) {gsub(/layer/,"top_layer")} \
                        else {gsub(/layer/,"middle_layer")}} \
          {print $0}' \
       $file  > $file.tmp \
	&& mv -f $file.tmp $file ) || echo "Failed translation in: $file"
	exit 0
fi
exit 1

