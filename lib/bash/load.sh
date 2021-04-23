sleep 1
#printf "LOADING..."
while [ 0 ]; do
    array_load=( "LOADING---" "LOADING.--" "LOADING..-" 'LOADING...' )
    for show_load in ${array_load[@]}; do
        #show_load="$str_load"
        number_of_char=`echo "$show_load" | wc -c`
        number_of_char=$(($number_of_char-1));
        printf "\b%.0s" `seq 1 $number_of_char`
        printf "$show_load";
        array_character=( "|" "/" "-" '.' )
        for i in ${array_character[@]}; do
            printf "%s" "${i}";
            sleep 0.2;
            printf "\b";
        done;
    done;
done;



# show_step="$i/$number_lines TE"
# number_of_char=`echo "$show_step" | wc -c`
# number_of_char=$(($number_of_char-1));
# printf "\b%.0s" `seq 1 $number_of_char`
# printf "$show_step";



