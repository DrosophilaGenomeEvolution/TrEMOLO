TE_POS_SIZE=$1
SAM_IN=$2
SAM_OUT=$3

awk -v size_windows_pos="15" -v max_size_TE_fk="30000"  -v size="10" 'BEGIN{
    qual=""
    nuc=""; 
    for(i=0; i<max_size_TE_fk; i++){
        nuc  = nuc"N"
        qual = qual"8"
    };

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
            if( caracter != "H" && caracter != "S" && caracter != "I" && number != ""){
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
            dico[$1, i] = $2":"$3
            i+=1
        }
        else{
            i=1
            chrom = $1
            dico[$1, i] = $2":"$3
            i+=1
        }
        next;
    }

OFS="\t"{ 
    cigar=$6
    size_seq=parse_cig_length(cigar)
    len=split($6, CIG,"S"); 
    n=match(CIG[1], /^[0-9]+$/); 

    #$3 chrom

    position_debut_read=(0 != 0)
    position_end_read=(0 != 0)

    i=1
    while( dico[$3, i] && ! position_debut_read && ! position_end_read ){
        split(dico[$3, i], pos_and_size, ":")
        position = pos_and_size[1]
        size_in  = pos_and_size[2]

        if((position - size_windows_pos <= $4 && $4 <= position + size_windows_pos)){
            position_debut_read=(position - size_windows_pos <= $4 && $4 <= position + size_windows_pos)
        }

        if((position - size_windows_pos <= $4 + size_seq && $4 + size_seq <= position + size_windows_pos)){
            position_end_read=(position - size_windows_pos <= $4 + size_seq && $4 + size_seq <= position + size_windows_pos)
        }
        #position_debut_read=(position - size_windows_pos <= $4 && $4 <= position + size_windows_pos)

        i += 1
    }

    # if( position_debut_read && position_end_read){
    #     print "kk", size_in, size_seq, position, $4, $6
    # }

    if( position_debut_read || position_end_read){
        
        if(len >= 2 && n != 0){
            rest=size_in-CIG[1]

            if( position_debut_read ){
                deb=size"M"(CIG[1]+rest)"I"; 
                if( size + rest < 0 ){
                    deb=((size+rest)*-1)"S"deb; 
                }
                $4=$4-size;
            }
            else {
                deb=CIG[1]"S"
            }

            match(CIG[2], /[A-Z][0-9]+$/);

            if(len == 3 && position_end_read){
                rest2 = size_in - substr(CIG[2], RSTART+1, RLENGTH)
                $6=deb""substr(CIG[2], 1, RSTART)""(substr(CIG[2], RSTART+1, RLENGTH)+rest2)"I"size"M";
                
                if( size + rest2 < 0){
                    $6=$6""( (size + rest2) * -1)"S"; 
                }
            }
            else{
                $6=deb""CIG[2]"S";
            }

            if($10 != "*"){
                if( size + rest > 0 && position_debut_read ){
                    seq=substr(nuc, 1, size + rest)
                    ql=substr(qual, 1, size + rest)
                    $10=seq""$10;
                    $11=ql""$11;
                }

                if( len == 3 && position_end_read ){
                    if(size + rest2 > 0){
                        seq=substr(nuc, 1, size + rest2)
                        ql=substr(qual, 1, size + rest2);
                        $10=$10""seq;
                        $11=$11""ql;
                    }
                }  
            }
        }
        else {
            #HARD
            len=split($6, CIG,"H"); 
            n=match(CIG[1], /^[0-9]+$/); 
            
            if(len >= 2 && n != 0){
                if( position_debut_read ){
                    deb=size"M"size_in"I";
                    $4=$4-size;
                }
                else{
                    deb=CIG[1]"H"
                }

                #$6=size"M"size_in"I"; 
                match(CIG[2], /[A-Z][0-9]+$/);


                if(len == 3 && position_end_read ){
                    $6=deb""substr(CIG[2], 1, RSTART)""size_in"I"size"M";
                }
                else {
                    $6=deb""CIG[2]"H";
                }

                
                if($10 != "*"){
                    seq=substr(nuc, 1, size + size_in)
                    ql=substr(qual, 1, size + size_in);
                    if( position_debut_read ){
                        $10=seq""$10;
                        $11=ql""$11;
                    }

                    if(len == 3  && position_end_read ){
                        $10=$10""seq;
                        $11=$11""ql;
                    }
                }
            }
        }
    }

    print $0;
}' ${TE_POS_SIZE} ${SAM_IN} \
> ${SAM_OUT}

