pid=`ps -f | grep "TrEMOLO/lib/bash/load.sh" | awk '$8=="bash" {print $2}'`
#echo "PID=${pid}"
kill -s SIGTERM ${pid} 
kill -s SIGINT ${pid}
echo 