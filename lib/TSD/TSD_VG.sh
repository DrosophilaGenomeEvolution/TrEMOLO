
INPUT_OUTPUT=$1

number_ok=`grep OK ${INPUT_OUTPUT} -c`
number_ko=`grep KO ${INPUT_OUTPUT} -c`
number_total=`grep ">" ${INPUT_OUTPUT} -c`
number_k_o=`grep "K-O" ${INPUT_OUTPUT} -c`

#RESUME
if [ -n "$number_total" ]; then
    
    echo "OK/total : $number_ok/$number_total" >> ${INPUT_OUTPUT}
    echo "KO/total : $number_ko/$number_total" >> ${INPUT_OUTPUT}
    echo "OK+KO/total : $(($number_ok+$number_ko))/$number_total" >> ${INPUT_OUTPUT}
    echo "K-O/total : $number_k_o/$number_total" >> ${INPUT_OUTPUT}
    echo "OK+K-O/total : $(($number_ok+$number_k_o))/$number_total" >> ${INPUT_OUTPUT}
    echo "OK% : $(($number_ok*100/$number_total))%" >> ${INPUT_OUTPUT}
else
    echo "ERROR : ${INPUT_OUTPUT}";
fi

