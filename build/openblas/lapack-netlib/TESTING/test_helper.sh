rm -f $3
$1 < $2
grep -q FATAL $3
if [ $? -eq 0 ]; then
echo Error
exit 1
else
exit 0
fi
