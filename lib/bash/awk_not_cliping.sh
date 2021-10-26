TE_POS_SIZE=$1
SAM_IN=$2
SAM_OUT=$3

awk -v size_windows_pos="30" -v max_size_TE_fk="30000"  -v size="100" 'BEGIN{
    qual=""
    nuc=""; 
    for(i=0; i<max_size_TE_fk; i++){
        nuc  = nuc"N"
        qual = qual"8"
    };
    print "" > "number_clipped.txt"
}

function parse_cig_length(cigar){
    sum=0
    number=""
    #print "len=" length(cigar)
    for(i=1; i<=length(cigar); i++){
        caracter=substr(cigar, i, 1)
        #print "c="caracter
        if( match(caracter, /^[0-9]$/) == 0 ){
            #print number""caracter
            if( caracter != "H" && caracter != "I" && caracter != "S" && number != ""){
                sum += number
            }
            number=""
        }
        else{
            number = number""caracter 
        }
    };
    return sum;
};

# Traitement le 1er fichier, on rempli un tableau
    
    NR == FNR {  
        if(chrom == $1){
            dico[$1, i] = $2":"$3":"$4
            i+=1
        }
        else{
            i=1
            chrom = $1
            dico[$1, i] = $2":"$3":"$4
            i+=1
        }
        next;
    }

OFS="\t"{ 
    cigar=$6
    size_seq=parse_cig_length(cigar)

    position_debut_read=(0 != 0)
    position_end_read=(0 != 0)

    size_in_deb = 0;
    size_in_end = 0;

    i=1
    while( dico[$3, i] && ( ! position_debut_read || ! position_end_read ) ){
        split(dico[$3, i], pos_and_size, ":")
        position = pos_and_size[1]
        size_in  = pos_and_size[2]
        TE_ID    = pos_and_size[3]

        if( ! position_debut_read && (position - size_windows_pos <= $4 && $4 <= position + size_windows_pos) && match($6, /^[0-9]+[SH]/) ){
            position_debut_read=(position - size_windows_pos <= $4 && $4 <= position + size_windows_pos)
            size_in_deb = pos_and_size[2]
        }

        if( ! position_end_read && (position - size_windows_pos <= $4 + size_seq && $4 + size_seq <= position + size_windows_pos)  && match($6, /[0-9]+[SH]$/) ){
            position_end_read = (position - size_windows_pos <= $4 + size_seq && $4 + size_seq <= position + size_windows_pos)
            size_in_end       = pos_and_size[2]
        }

        i += 1
    }

    print $1":DEBU_CIG:"$6 >> "number_clipped.txt";

    if( position_debut_read || position_end_read){
        print $1":"$3":"position":"size":"TE_ID":SOFT:"size_in_deb >> "number_clipped.txt"

        if( match($6, /[0-9]+S/) ){
            rest = 0;
            
            #DEB
            
            if( position_debut_read && match($6, /^[0-9]+S/) ){
                SOFT_DEB = substr($6, RSTART, RLENGTH-1)
                rest     = size_in_deb - SOFT_DEB

                print "--position_debut_read:"position_debut_read"; rest:"rest" ; SOFT_DEB:"SOFT_DEB >> "number_clipped.txt"

                deb=size"M"( SOFT_DEB + rest )"I";
                 
                if( size + rest < 0 ){
                    deb=((size + rest)*-1)"S"deb;
                }
                $4=$4 - size;
            }
            else if( match($6, /^[0-9]+S/) ){
                deb=substr($6, RSTART, RLENGTH)
            }
            else{
                deb=""
            }


            rest2 = 0;

            #END
            print "position_end_read:"position_end_read"; rest2:"rest2" ; SOFT_END:"SOFT_END"; size_in_end:"size_in_end >> "number_clipped.txt"

            if( match($6, /[0-9]+S$/) && position_end_read ){
                SOFT_END = substr($6, RSTART, RLENGTH-1)
                rest2    = size_in_end - SOFT_END

                print "::!!position_end_read:"position_end_read"; rest2:"rest2" ; SOFT_END:"SOFT_END"; size_in_end:"size_in_end >> "number_clipped.txt"

                match($6, /[0-9A-RT-Z]+[A-RT-Z]/);

                $6=deb""substr($6, RSTART, RLENGTH)""(SOFT_END + rest2)"I"size"M";

                if( size + rest2 < 0 ){
                    $6=$6""( (size + rest2) * -1 )"S"; 
                }
            }
            else{
                
                if( ! position_debut_read ){
                    match($6, /[0-9A-Z]+/);
                    $6=deb""substr($6, RSTART, RLENGTH)
                }
                else{
                    match($6, /[0-9A-RT-Z]+[0-9S]*$/);
                    $6=deb""substr($6, RSTART, RLENGTH)
                }

            }


            if($10 != "*"){
                print "10:OK; size:"size"; " >> "number_clipped.txt"
                if( size + rest > 0 && position_debut_read ){
                    print "DEB_ADD" >> "number_clipped.txt"
                    seq=substr(nuc, 1, size + rest)
                    ql=substr(qual, 1, size + rest)
                    $10=seq""$10;
                    $11=ql""$11;
                }

                if( position_end_read ){
                    if( size + rest2 > 0 ){
                        print "END_ADD" >> "number_clipped.txt"
                        seq=substr(nuc, 1, size + rest2)
                        ql=substr(qual, 1, size + rest2);
                        $10=$10""seq;
                        $11=$11""ql;
                    }
                }  
            }
        }

        print "CIG_SOFT:"$6 >> "number_clipped.txt";
        
        ##
        if( match($6, /[0-9]+H/) ) {
            print $1":"$3":"position":"size":"TE_ID":HARD" >> "number_clipped.txt"
            #HARD
            
            if( match($6, /^[0-9]+H/) ){
                HARD_DEB = substr($6, RSTART, RLENGTH-1)

                print "--position_debut_read:"position_debut_read"; HARD_DEB:"HARD_DEB >> "number_clipped.txt"

                if( position_debut_read ){
                    deb=size"M"size_in_deb"I";
                    $4=$4-size;
                }
                else{
                    deb=HARD_DEB"H"
                }
            }
            else{
                deb=""
            }

            #END
            if( match($6, /[0-9]+H$/) ){
                HARD_END = substr($6, RSTART, RLENGTH-1)

                print "!!position_end_read:"position_end_read"; HARD_END:"HARD_END >> "number_clipped.txt"

                match($6, /[0-9A-GI-Z]+[A-GI-Z]/);
                if( position_end_read ){
                    $6=deb""substr($6, RSTART, RLENGTH)""size_in_end"I"size"M";
                }
                else {
                    $6=deb""substr($6, RSTART, RLENGTH)""HARD_END"H";
                }
            }
            else{
                HARD_END = substr($6, RSTART, RLENGTH-1)
                match($6, /[0-9A-GI-Z]+[A-GI-Z]/);
                if( HARD_END != "" ){
                    $6=deb""substr($6, RSTART, RLENGTH)""HARD_END"H";
                }
                else{
                    $6=deb""substr($6, RSTART, RLENGTH);
                }
            }

            print "CIG_HARD:"$6 >> "number_clipped.txt";

            if($10 != "*"){
                
                if( position_debut_read ){
                    seq=substr(nuc, 1, size + size_in_deb)
                    ql=substr(qual, 1, size + size_in_deb);
                    $10=seq""$10;
                    $11=ql""$11;
                }

                if( position_end_read ){
                    seq=substr(nuc, 1, size + size_in_end)
                    ql=substr(qual, 1, size + size_in_end);
                    $10=$10""seq;
                    $11=$11""ql;
                }
            }
        }
    }

    print $0;
} 

END{
    print "FINISH" >> "number_clipped.txt"
}' ${TE_POS_SIZE} ${SAM_IN} \
> ${SAM_OUT}

