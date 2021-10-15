# Practical Attack

The following scripts/sources were used to collect the data for the evaluation of practical attacks. The most important files are:

* ```jupyter/sefa.ipynb```: contains the python code that was used to communicate and tamper with the chipwhisperer victim device.
* ```simpleserial-sefa/simpleserial-sefa.c```: contains the main firmware code that was used for evaluating the attacks.

We also provide compiled binaries and the makefiles that we have used. In order to use this code simply setup an ordinary installation of the CW framework and copy the provided files to their usual locations. To reprocue the fault attacks it will very likely be necessary to adjust the parameters of the clock glitch as they will differ slightly between each board.
Also, ```simpleserial-sefa/simpleserial-sefa.c``` contains the implementations of both the SIFA-protected variant of the Keccak S-box and the ordinary DOM variant. Simply comment in the one that you want to evaluate.

