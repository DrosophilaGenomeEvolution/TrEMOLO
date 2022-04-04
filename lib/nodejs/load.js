
const tty = require('tty');
const status = process.stdout.isTTY;

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

strand     = 1;
rand_dep = getRandomInt(7) + 1;
rand_end = getRandomInt(7);
GENOME   = "________    ________";

//❛ั∗✺✖

async function load() {

    distance = Math.abs(rand_dep - (rand_end + 14));
    e=0;
    while (1) {
        if (show) {
            
            if ( e >= distance ) {
                e=0;
                strand *= -1;
                rand_dep = 20 - distance - rand_dep;
                rand_end = getRandomInt(7);
                distance = Math.abs(rand_dep - (rand_end + 15));
            }
            process.stdout.cursorTo(0, `${process.stdout.rows}`-2)
            process.stdout.clearLine(1)
            
            TE = GENOME.replace(/./g, " ");
            TE = TE.replaceAt(rand_dep + e - 1, "∗✺∗")

            //
            if ( e == distance - 1 || e == 0 ) {
                TE = "";
            }

            //STRAND TE MOVEMENT
            if ( strand == 1 ) {
                process.stdout.write("    " + TE);
            }
            else{
                process.stdout.write("    " + TE.split("").reverse().join(""));
            }
            
            process.stdout.cursorTo(0, `${process.stdout.rows}`-1)

            //GENOME
            if ( e == distance - 1 ){
                if (strand == 1) {
                    genome = GENOME.replaceAt(rand_dep + e - 1, "✺")
                }
                else{
                    genome = GENOME.replaceAt(rand_dep + e - 1, "✺").split("").reverse().join("")
                }
                process.stdout.write("    " + genome + "   ");
            }
            else if( e == 0 ){
                if (strand == 1) {
                    genome = GENOME.replaceAt(rand_dep + e, "∗✺∗")
                }
                else{
                    genome = GENOME.replaceAt(rand_dep + e, "∗✺∗").split("").reverse().join("")
                }
                process.stdout.write("    " + genome + "   ");
            }
            else{
                process.stdout.write("    " + GENOME + "   ");
            }

            process.stdout.cursorTo(0)
            if ( e == distance - 1 || e == 0 ){
                await sleep(300);
            }
            else{
                await sleep(35);
            }
            
            process.stdout.clearLine(-1)
            process.stdout.cursorTo(0)
            e+=1;
        }
        else{
            await sleep(1000);
        }
    }
}

//console.log(tty.isatty(1))

if ( status ) {
    /*process.stdout.on('resize', () => {
        console.log('screen size has changed!');
        console.log(`${process.stdout.columns}x${process.stdout.rows}`);
    });*/

    //process.stdin.setRawMode(true);

    console.log("MODE TERMINAL");

    //SIGINT   2
    process.on('SIGINT', () => {
        process.stdout.cursorTo(0, `${process.stdout.rows}`-3)
        process.stdout.clearLine(1)
        process.stdout.cursorTo(0, `${process.stdout.rows}`-2)
        process.stdout.clearLine(1)
        process.stdout.cursorTo(0, `${process.stdout.rows}`-1)
        process.stdout.clearLine(1)
        console.log("WAITING DONE !");
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
            console.log(" ••• WAIT TE ARE MOVING •••\n");
            show = true;
        }
    });

    process.on('SIGTERM', () => {
        process.stdout.clearLine(-1)
        process.stdout.cursorTo(0)
        console.log("WAITING DONE !");
        process.exit(0);
    });

    process.on('SIGQUIT', () => {
        process.stdout.clearLine(-1)
        process.stdout.cursorTo(0)
        console.log("Quit");
        process.exit(0);
    });
    
    load();
}
else {
    console.log("REDIRECTING MODE");
}

