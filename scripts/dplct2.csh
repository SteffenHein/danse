#!/bin/csh
set init=1
set final=400

set COPIED_FILE="_TU"
set TARGET_FILE="_TU"
set SUFFIX=""
# now iterate

@ i = $init
while ( $i <= $final )
	if (-f $i$COPIED_FILE ) then
		echo "copying file $i$COPIED_FILE to ./tmp/$i$TARGET_FILE"
		cp $i$COPIED_FILE "./tmp/"$i$TARGET_FILE
	endif
	@ i ++
end
#
echo "Done"
