
const tty = require('tty');
const status = process.stdout.isTTY;

//const readlineModule = require('readline');

var show = false;

function getRandomInt(max) {
    return Math.floor(Math.random() * max);
}

function sleep(ms) {
    return new Promise((resolve) => {
        setTimeout(resolve, ms);
    });
}

String.prototype.replaceAt = function(index, replacement) {
    return this.substring(0, index) + replacement + this.substring(index + replacement.length);
}

function init_TE(lvl = 0){
    let TE = {
        pat: 0,
        strand: 1,
        distance: 0,
        rand_dep: getRandomInt(SIZE_CHROM) + 1,
        rand_end: getRandomInt(SIZE_CHROM),
        level: lvl
    }
    return TE;
}

GENOME   = "_______________    _______________";

SIZE_GENOME_EMPTY = GENOME.match(/^_+\s+/)[0].length - 1;
SIZE_CHROM        = GENOME.match(/^_+/)[0].length - 1;

genome_style = ["-", "•", "_", "✖"];
num_style    = getRandomInt(genome_style.length);
GENOME       = GENOME.replace(/_/g, genome_style[num_style]);
//❛ั∗✺✖

var index_TE = 0;
var TEs      = [];
TEs.push(init_TE());

//console.log(TEs);

async function load() {

    TEs[index_TE].distance = Math.abs(TEs[index_TE].rand_dep - (TEs[index_TE].rand_end + SIZE_GENOME_EMPTY));
    TEs[index_TE].pat = 0;
    while (1) {
        if (show) {
            for( index_TE = 0; index_TE < TEs.length; index_TE++ ){
                if ( TEs[index_TE].pat >= TEs[index_TE].distance ) {
                    TEs[index_TE].pat=0;
                    TEs[index_TE].strand  *= -1;
                    TEs[index_TE].rand_dep = GENOME.length - TEs[index_TE].distance - TEs[index_TE].rand_dep;
                    TEs[index_TE].rand_end = getRandomInt(SIZE_CHROM);
                    TEs[index_TE].distance = Math.abs(TEs[index_TE].rand_dep - (TEs[index_TE].rand_end + SIZE_GENOME_EMPTY + 2));
                }
                process.stdout.cursorTo(0, `${process.stdout.rows}`-2 - TEs[index_TE].level)
                process.stdout.clearLine(1)
                
                TE = GENOME.replace(/./g, " ");
                TE = TE.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat - 1, "∗✺∗")

                //
                if ( TEs[index_TE].pat == TEs[index_TE].distance - 1 || TEs[index_TE].pat == 0 ) {
                    TE = "";
                }

                //STRAND TE MOVEMENT
                if ( TEs[index_TE].strand == 1 ) {
                    process.stdout.write("    " + TE);
                }
                else{
                    process.stdout.write("    " + TE.split("").reverse().join(""));
                }
                
                process.stdout.cursorTo(0, `${process.stdout.rows}`-1)

                //GENOME
                if ( TEs[index_TE].pat == TEs[index_TE].distance - 1 ) {
                    if ( TEs[index_TE].strand == 1 ) {
                        genome = GENOME.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat - 1, "✺")
                    }
                    else {
                        genome = GENOME.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat - 1, "✺").split("").reverse().join("")
                    }
                    process.stdout.write("    " + genome + "    ");
                }
                else if ( TEs[index_TE].pat == 0 ) {
                    if (TEs[index_TE].strand == 1) {
                        if( TEs[index_TE].rand_dep >= 0 )
                            genome = GENOME.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat, "∗✺∗")
                        else
                            genome = GENOME.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat + 1, "∗✺∗")
                    }
                    else {
                        genome = GENOME.replaceAt(TEs[index_TE].rand_dep + TEs[index_TE].pat, "∗✺∗").split("").reverse().join("")
                    }
                    process.stdout.write("    " + genome + "    ");
                }
                else {
                    process.stdout.write("    " + GENOME + "    ");
                }

                process.stdout.cursorTo(0)
                if ( TEs[index_TE].pat == TEs[index_TE].distance - 1 || TEs[index_TE].pat == 0 ) {
                    await sleep(300 - ( TEs.length  * 10 ) );
                }
                else {
                    await sleep(35 - ( TEs.length  * 5 ) );
                }
                
                process.stdout.clearLine(-1);
                process.stdout.cursorTo(0);
                TEs[index_TE].pat += 1;
            }
        }
        else{
            await sleep(1000);
        }
    }
}

//console.log(tty.isatty(1))

if ( status ) {

    //readlineModule.emitKeypressEvents(process.stdin);
    //process.stdin.setRawMode(true);

    //RESIZE SCREEN
    process.stdout.on('resize', () => {

        show = false;
        //console.log('screen size has changed!');
        //console.log(`${process.stdout.columns}x${process.stdout.rows}`);
        GENOME = ""
        for (var i = 0; i < ((`${process.stdout.columns}` - 9)/2)-2; i++) {
            GENOME += "_";
        }

        GENOME += "    ";

        for (var i = 0; i < ((`${process.stdout.columns}` - 9)/2)-2; i++) {
            GENOME += "_";
        }

        SIZE_GENOME_EMPTY = GENOME.match(/^_+\s+/)[0].length - 1;
        SIZE_CHROM        = GENOME.match(/^_+/)[0].length - 1;

        GENOME            = GENOME.replace(/_/g, genome_style[num_style])
        
        for(var tmp_index_TE = 0; tmp_index_TE < TEs.length; tmp_index_TE++ ){
            TEs[tmp_index_TE].rand_dep = getRandomInt(SIZE_CHROM) + 1;
            TEs[tmp_index_TE].rand_end = getRandomInt(SIZE_CHROM);
            TEs[tmp_index_TE].distance = Math.abs(TEs[tmp_index_TE].rand_dep - (TEs[tmp_index_TE].rand_end + SIZE_GENOME_EMPTY + 2));
            TEs[tmp_index_TE].pat = 0;
        }

        show = true;
    });

    //Type Enter
    process.stdin.on('keypress', function (letter, key) {
        //console.log(letter, key)
        if ( key && key.name == "return" && TEs.length < 5 ) {
            console.log("");
            TEs.push(init_TE(TEs.length));
        }
        else if ( key.ctrl && ( key.name == "c" || key.name == "C" ) ) {
            process.exit(0);
        }
    });

    //process.stdin.setRawMode(true);

    console.log("LOADING ACTIVATE");

    //SIGINT   2
    process.on('SIGINT', () => {
        process.stdout.cursorTo(0, `${process.stdout.rows}`-3)
        process.stdout.clearLine(1)
        process.stdout.cursorTo(0, `${process.stdout.rows}`-2)
        process.stdout.clearLine(1)
        process.stdout.cursorTo(0, `${process.stdout.rows}`-1)
        process.stdout.clearLine(1)
        console.log("ANOMALY DETECTED");
        process.exit(0);
    });

    //SIGUSR2   12    
    process.on('SIGUSR2', () => {
        if ( show ) {
            show = false;
            process.stdout.cursorTo(0, `${process.stdout.rows}`-3)
            process.stdout.clearLine(1)
            process.stdout.cursorTo(0, `${process.stdout.rows}`-2)
            process.stdout.clearLine(1)
            process.stdout.cursorTo(0, `${process.stdout.rows}`-1)
            process.stdout.clearLine(1)
            console.log("WAITING DONE !");
        }
    });

    //SIGUSR1   10
    process.on('SIGUSR1', () => {
        if ( ! show ) {
            process.stdout.cursorTo(0)
            process.stdout.clearLine(1)
            console.log("")
            console.log(" ••• WAIT TE ARE MOVING •••\n");
            show = true;
        }
    });

    process.on('SIGTERM', () => {
        process.stdout.clearLine(-1)
        process.stdout.cursorTo(0)
        console.log("BYE !!");
        process.exit(0);
    });

    process.on('SIGQUIT', () => {
        process.stdout.clearLine(-1)
        process.stdout.cursorTo(0)
        console.log("Quit");
        process.exit(0);
    });
    
    console.log("")
    load();
}
else {
    console.log("LOADING DEACTIVATE");
}

