trap 'printf "\b%.0s" `seq 1 46` >&2; echo -e "\nEND LOAD\n" >&2; exit 1' SIGTERM SIGINT
sleep 1
#printf "LOADING..."
chain="...................................                            \b\b\b\b\b\b\b\b\b\b ";
while [ 0 ]; do
    array_load=( 'LOADING----' 'LOADING.---' 'LOADING..--' 'LOADING...-' 'LOADING....' )
    for show_load in ${array_load[@]}; do
        #show_load="$str_load"
        number_of_char=`echo "$show_load" | wc -c`
        number_of_char=$(($number_of_char-1));
        printf "\033[s"
        #printf "\033[%d;%dH\b" "1" "1";
        
        if [ `echo $chain | grep -o ".$"` = ">" ]; then
            chain=`echo -e $chain | sed '0,/>/s//./'`
        else
            chain=`echo -e $chain | sed '0,/\./s//>/'`
        fi
        

        #printf "\b%.0s" `seq 1 $number_of_char`
        printf "%s" "$show_load $chain" ;
        number_of_char=`echo "$show_load" | wc -c`
        number_of_char=$(($number_of_char-1));
        printf "\033[u";
        sleep 0.2
        # printf "\033[s"
        # printf "\033[%d;%dH\b" "1" "1";
        # printf "\b%.0s" `seq 1 $number_of_char`
        # printf "\033[K";
        # printf "\033[u";

    done;
done >&2;



# show_step="$i/$number_lines TE"
# number_of_char=`echo "$show_step" | wc -c`
# number_of_char=$(($number_of_char-1));
# printf "\b%.0s" `seq 1 $number_of_char`
# printf "$show_step";



# printf "\033[39m\033[39m - Reset colour"
# printf "\033[2K - Clear Line"
# printf "\033[<L>;<C>H OR \033[<L>;<C>f puts the cursor at line L and column C."
# printf "\033[<N>A Move the cursor up N lines"
# printf "\033[<N>B Move the cursor down N lines"
# printf "\033[<N>C Move the cursor forward N columns"
# printf "\033[<N>D Move the cursor backward N columns"
# printf "\033[2J Clear the screen, move to (0,0)"
# printf "\033[K Erase to end of line"
# printf "\033[s Save cursor position"
# printf "\033[u Restore cursor position"
# printf " "
# printf "\033[4m  Underline on"
# printf "\033[24m Underline off"
# printf "\033[1m  Bold on"
# printf "\033[21m Bold off"

