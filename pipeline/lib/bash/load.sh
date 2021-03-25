sleep 1
printf "LOADING..."
while [ 0 ]; do
    array_character=( "|" "/" "-" '.' )
    for i in ${array_character[@]}; do
        printf "%s" "${i}";
        sleep 0.2;
        printf "\b";
    done;
done;
