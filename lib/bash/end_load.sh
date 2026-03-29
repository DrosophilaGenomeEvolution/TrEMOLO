#!/usr/bin/env bash

set -E

trap 'status=$?; printf "%s\n" "[ERROR][end_load] line=$LINENO cmd=$BASH_COMMAND" >&2; exit $status' ERR

#pid=`ps -f | grep "TrEMOLO/lib/bash/load.sh" | awk '$8=="bash" {print $2}'`
pid=`ps -f | grep "TrEMOLO/lib/nodejs/load.js" | awk '$8=="node" {print $2}'`
if [ -n "${pid}" ]; then
    kill -s 12 ${pid} 2>/dev/null || true
fi;
