#!/bin/csh
set init=0
set final=500

set COPIED_FILE="flow1"
set TARGET_FILE="10.E-4_seconds:_"
set SUFFIX=""
# now iterate

@ i = $init
while ( $i <= $final )
	if (-f ./$COPIED_FILE$i) then
		echo "copying file $COPIED_FILE$i to $TARGET_FILE$i$SUFFIX"
		cp $COPIED_FILE$i $TARGET_FILE$i$SUFFIX
	endif
	@ i ++
end
#
echo "Done"
