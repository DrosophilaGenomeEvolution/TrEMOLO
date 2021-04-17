get_seq_with_id () {
    grep "^[^>][A-Za-z\.]*[0-9]*[0-9]" $1 -B 1 | grep -v "\-\-" | grep "[0-9]" | awk 'BEGIN{head=""}
    {
        if( substr($0, 1, 1) != ">" ){

            #>Assemblytics_b_1220:+:Repeat_expansion::3L_RaGOO_RaGOO_RaGOO:10739030-10763902
            size_sp = split(head, head_sp, ":")
            split(head_sp[6], coord, "-")

            match($0, /[0-9]+[A-Z]+/);
            seq[0]  = substr($0, RSTART, RLENGTH)

            start = coord[1] + RSTART - 2 #-2 pour bedtools
            end   = coord[1] + RSTART + RLENGTH
            
            position_next = RSTART + RLENGTH
            

            match(seq[0], /[0-9]+/)
            ID = substr(seq[0], RSTART, RLENGTH)

            head=head_sp[1]":"head_sp[2]":"head_sp[3]":"head_sp[4]":"head_sp[5]":"start"-"end":"ID
            print head"\n"seq[0]

            
            n = match(substr($0, position_next, length($0)), /[0-9]+[A-Z]+/); 

            i = 1; 
            while( substr($0, position_next + RSTART, RLENGTH) != substr($0, position_next + RSTART, length($0)) && n != "" ){
                seq[i]=substr(substr($0, position_next, length($0)), RSTART, RLENGTH)
                
                start = coord[1] + position_next - 2 #-1 pour bedtools
                end   = coord[1] + position_next + RLENGTH
                
                match($0, seq[i]); 
                position_next = RSTART + RLENGTH

                match(seq[i], /[0-9]+/)
                ID   = substr(seq[i], RSTART, RLENGTH)

                head = head_sp[1]":"head_sp[2]":"head_sp[3]":"head_sp[4]":"head_sp[5]":"start"-"end":"ID
                print head"\n"seq[i]
                
                n = match(substr($0, position_next, length($0)), /[0-9]+[A-Z]+/); 

                i += 1
                
            } 

            if( n != 0 ){
                start = coord[1] + position_next - 2
                end   = coord[1] + position_next + RLENGTH

                match(substr($0, position_next, length($0)), /[0-9]+/)
                ID = substr(substr($0, position_next, length($0)), RSTART, RLENGTH)

                head=head_sp[1]":"head_sp[2]":"head_sp[3]":"head_sp[4]":"head_sp[5]":"start"-"end":"ID
                print head"\n"substr($0, position_next, length($0))
            }
        }
        else{
            head=$0
        }

    }' | grep "[0-9]" > $2
}

echo "GET"
get_seq_with_id $1 $2