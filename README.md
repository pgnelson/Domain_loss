# Domain_loss
This repository contains the code used in the domain loss rate over evolutionary time project

These scripts were written to allow high throughput processing using a pbs system of non-coding sequences to identify domains that may be missing from current annotation.

"run_aa_processes.py" initiated and controlls the process. It calls "generate_aa_files.py" to generate amino acid files, and then several threads of "process_aa_files.py". It also cleans up directories that have been processed and collates results.
If there are no amino acid directories present it works best to call "generate_aa_files" first, wait for at least 100 directories of amino acid files to be generated, and then call "run_aa_processes".
