# Modified EQL/EVP
Built from EQL/EVP (Risacher and Clement 2001), a program for simulation of evaporation of waters to high salinities. The goal of this project is to turn EQLEVP into functions that can be scripted, rather than an interactive prompt interface, in order to run iterative simulations. The original FORTRAN90 programs run as interactive, command-line prompts that ask the user for input before proceeding; this makes iterative calls to the program a slow process that requires user attention.

Currently, the original FORTRAN90 code has been translated to Python and made into an object-oriented/functional structure that can be easily scripted for iterative and batch runs. The program runs rather slowly in Python, so I am working on a C++ translation in addition to the Python version. Again, this will be object-oriented, which will make for easy and quick simulations.
