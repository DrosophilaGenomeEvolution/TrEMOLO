trap 'printf "\b%.0s" `seq 1 46` >&1; echo -e "\nEND LOAD\n" >&1; exit 1' SIGTERM SIGINT;
sleep 1;
chain="•••";
#echo ".·’·.•’•.٭*"
while [ 0 ]; do
    array_load=( 'LOADING.    ' 'LOADING •   ' 'LOADING  *  ' 'LOADING   • ' 'LOADING    .' 'LOADING   • ' 'LOADING  *  ' 'LOADING •   ' 'LOADING.    ')
    for ((i = 0; i < ${#array_load[@]}; i++));
    do
        show_load="${array_load[$i]}"

        number_of_char=`echo "$show_load" | wc -c`
        number_of_char=$(($number_of_char-1));
        printf "\033[s";
        
        if [ `echo $chain | grep -e "[.•]*" -o | grep -o ".$"` = "•" ]; then
            chain=`echo -e "$chain" | sed '0,/•/s//./'`
        else
            chain=`echo -e "$chain" | sed '0,/\./s//•/'`
        fi;
        
        printf "%s" "$show_load $chain" ;
        number_of_char=`echo "$show_load" | wc -c`
        number_of_char=$(($number_of_char-1));
        printf "\033[u";
        sleep 0.2;
    done;
done >&1;
