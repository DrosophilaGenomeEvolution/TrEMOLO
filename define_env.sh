#Define environement shell utils

#COLOR#
##############
CYAN='\033[96m'
VIOLET='\033[95m'
RED='\033[91m'
BLUE='\033[94m'
CYAN='\033[96m'
GREEN='\033[92m'
WARNING='\033[93m'
FAIL='\033[91m'
BOLD='\033[1m'
END='\033[0m'
UNDERLINE='\033[4m'


#FUNCTION#
#################

run_cmd () { 
    cmd=$1;
    log=$2;
    title=$3;
    show_out=$4;

    #printf "\n%s\n\n" "${CYAN} [SNK]--[`date`] $title ${END}"
    printf "%s\n" "$cmd";
    (eval "$cmd" 2>> ${log}.err 1>> ${log}.out && echo -e "${GREEN} TASK $title is DONE ${END} check ${log}.out AND ${log}.err") \
        || echo -e "${RED} TASK ERROR ${END} :  $cmd  \n PLEASE CHECK : ${log}.err";
    
    if [ $show_out = "True" ]; then
        cat ${log}.out;
    fi;
};


message () {
    message=$1

    #echo -e "${CYAN} ${message} ${END}"
    printf "%b" "${CYAN} ${message} ${END}"
}


message_fail () {
    message=$1

    #echo -e "${RED} ${message} ${END}"
    printf "%b" "${RED} ${message} ${END}"
}


begin_load (){
    (bash ${path_to_pipline}/lib/bash/load.sh &) || echo
}


end_load (){
    sh ${path_to_pipline}/lib/bash/end_load.sh;
}
