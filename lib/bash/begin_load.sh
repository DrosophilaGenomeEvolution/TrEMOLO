pid=`ps -f | grep "TrEMOLO/lib/nodejs/load.js" | awk '$8=="node" {print $2}'`
if [ -n "${pid}" ]; then
    kill -s 10 ${pid} 
fi;