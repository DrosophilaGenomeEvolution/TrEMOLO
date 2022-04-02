
const tty = require('tty');
const status = process.stdout.isTTY;

var show = false;

function sleep(ms) {
    return new Promise((resolve) => {
        setTimeout(resolve, ms);
    });
}

async function load() {
    //array_TE=['LOADING.    ', 'LOADING •   ', 'LOADING  *  ', 'LOADING   • ', 'LOADING    .', 'LOADING   • ', 'LOADING  *  ', 'LOADING •   ', 'LOADING.    '];
    array_TE     = ['_______      _______', '_______      _______', '___•___      _______', '____•__      _______', '______.      _______', '_______  .   _______', '_______     ._______', '_______      _•_____', '_______      ___•___', '_______      _______', '_______      _______'];
    array_GENOME = ['_______      _______', '___•___      _______', '_______      _______', '_______      _______', '_______      _______', '_______      _______', '_______      _______', '_______      _______', '_______      _______', '_______      ___•___', '_______      _______'];
    e=0;
    while (1) {
        if (show) {
            if ( e >= array_TE.length ) {
                e=0;
                array_TE.reverse();
                array_GENOME.reverse();
            }
            process.stdout.cursorTo(0, `${process.stdout.rows}`-2)
            process.stdout.clearLine(1)
            process.stdout.write("    " + array_TE[e].replace(/_/g, " "));
            process.stdout.cursorTo(0, `${process.stdout.rows}`-1)
            process.stdout.write("    " + array_GENOME[e]);
            process.stdout.cursorTo(0)
            await sleep(100);
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
        //process.stdout.clearLine(-1)
        //process.stdout.cursorTo(0)
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

