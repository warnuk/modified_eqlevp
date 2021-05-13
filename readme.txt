EQL/EVP programs
version 99/01
date: 03/1999




INSTALLATION AND USE
====================

WINDOWS 95/98, NT

- unzip the archive evp.zip to a specific directory

- to run the program, 2 options:

	+ double-click on the executable eql.exe

	+ or open a command window
		 go to the directory containing EVP files
		 type "eql"

remarks: It's a good idea to change the properties of the command window:
			- in Windows 95/98: set 50 lines to the command window 
			- in windows NT: increase the number of buffered lines. 
			  This will allow the user to scroll up and down in the window
			  and to see the outputs from the program.
			  
		 Do not forget that the Test data are in the test directory


UNIX

- uncompress the archive evp.zip (e.g. unzip evp.zip) and store the files 
in a specific directory (the main directory)

- replace the source files main_evp.f90 and main_eql.f90 by
main_evp_unix.f90 and main_eql_unix.f90

- compile both sources
  (e. g. under Solaris:  f90 -fast main_eql.f90 -o eql.exe
                         f90 -fast main_eqv.f90 -o eqv.exe  )

- run the program by typing "eql"

remark: Do not forget that the Test data are in the test directory



DESCRIPTION OF SOME INPUT PARAMETERS
====================================

Dilution factor:
				Especially useful if the initial solution is 
				oversaturated with respect to several minerals
				constraining the activity of water.
				The user can dilute the solution by giving the factor of
				dilution (format: 1.xx where xx represents the percent
				of dilution). The process is repeated until EQL receives
				a carriage return as an answer for the dilution factor.
				
Minerals data base modifications:
				The user can introduce or remove one or more minerals
				from the current set of minerals used in the simulation.
				
"Open" and "closed" system:
				In a "closed" system (equilibrium mode), minerals are 
				allowed to	redissolve. In an "open" system (fractional 
				crystallization mode, the amounts of precipitated
				minerals are removed at each step from the solution.
				
Automatic increment:	
				The increment is the percentage of water removed at each
				calculation step. Typically, values between 5% in dilute
				solutions and 0.1% in brines may be used in most 
				calculations.

Increment			
				The option "automatic" lets the program choose a 
				variable increment as a function of salinity. If the 
				option "manual" is on, the increment entered by the user 
				is constant throughout the simulation.		
				
End point:
				The user can fix a final salinity (in mg/l) for the
				simulation.
				
Screen output step:
				The progress of the evaporation is printed on the
				screen. By default, every increment is printed out
				(step=1). The user can choose an other step.
				
Output screen logging:
				If necessary, the data displayed on screen can also be 
				written in a file (its name is sample_name.log)
				
Results stored into files:
				3 files are created on demand. The user can define the
				name of the files or accept the standard names. 
				An increment step must be specified for the storage 
				(5 or 10 are good values). 
				
				The first file (default=sample_name.j[oc]&) contains
				for each saved increment the set of precipitated
				minerals and the current value of the parameters of the
				solution (fc: concentration factor, eva: amount of 
				evaporated water, ds: density, ph, alk, the 					
				concentrations of aqueous species, and total salinity). 			
				In addition, the program stores the composition of 
				the solution each time a paticular event occurs 
				during the simulation. Units are millimoles/l or
				millimoles/kg(H2O)
				
				The second file (default=sample_name.j[oc]@) contains
				the events that occurred during the simulation (start
				and end of dissolution or precipitation of minerals, 
				invariant points).
				
				The third file (default=sample_name.j[oc]%) contains
				for each saved increment the amount of precipitated
				minerals in moles.

Limits of convergence: 
				The Newton-Raphson method is an iterative calculation 
				that ends when all elements of a vector are lower than 
				an arbitrary small value: the limit of convergence. When 
				divergence occurs, one may try to lower the limit and 
				restart the simulation.
				


CONTENT OF evp.zip
==================

- README.TXT (this file)

- Program files (for Windows)

EQL.EXE: equilibrium program (to be used first)

EVP.EXE: evaporation program (launched by EQL.EXE)


- Data files

AQU.DAT : Names, atomic weights and electrical charges of the species
considered in EQL

AQUV.DAT : Names, atomic weights and electrical charges of the species
considered in EVP

COEFFT4 : Coefficients of  fourth degree polynomials for calculating
Pitzer's and silica parameters.

COMPLEX3 : Coefficients of fourth degree polynomials for calculating
dissociation constants of water, carbonate species, ion pairs and
complexes.

DENSITE : Coefficients for calculating the approximate density of the
solution.

MATRICE1 : Coefficients of partial derivatives of the EQP equation set
to be solved by the Newton-Raphson method.

MATRICE2 : Coefficients of partial derivatives of the EVP equation set
to be solved by the Newton-Raphson method.

MURTF2 : EQL mineral data base

MURTF3 : EVP mineral data base

MURTF0 : Temporary file created if the EVP mineral data base has been
modified.


- Test files (ASCII data files)

TEST_1.DAT
TEST_2.DAT

Default values of parameters used for the tests (unless otherwise specified):

	-input by file (f)
	-current directory (i.e. the directory containing the data files and 
		programs)
	-no change for temperature and log pCO2
	-initial solution not printed
	-no dilution
	-minerals data base not changed
	-closed system (c)
	-automatic increment 
	-no end point
	-screen output step 10
	-no output screen logging
	-results of the simulation stored in files
	-storage step 10
	-units = molarity
	-standard file names accepted
	-No modifications for the limits of convergence

How to run the tests:
	
	At the question: Input file/keyboard (f/k):
	answer:  f  or  press <enter>.

	At the question: Directory (current directory = <enter>):
	answer:  test

	At the question: File Name:
	answer:  test_1.dat   or   test_2.dat

	At the question: Number or name of analysis
	answer:  the name or the preceding mumber displayed on screen

In TEST_1.DAT

	Sample 1: CAL-19
	CAUTION: Log PCO2 = -3.5
	Results in	cal-19.jc%
			cal-19.jc&
			cal-19.jc@

	
	Sample 2: SAL-7
	CAUTION: Factor of dilution = 1.02
	Results in	sal-7.jc%
			sal-7.jc&
			sal-7.jc@


	Sample 3: UA-600
	CAUTION: 	increment = 0.1 
			storage step = 5
	Results in	ua-600.jc%
			ua-600.jc&
			ua-600.jc@
			

	Sample 4: UA-600
	CAUTION: Open system
	Results in	ua-600.jo%
			ua-600.jo&
			ua-600.jo@

	Sample 5: ASC-14
	Results in	asc-14.jc%
			asc-14.jc&
			asc-14.jc@


In TEST_2.DAT

	Sample 1: KMG_EUT
	CAUTION: the results of the simulation are stored in files with
		-storage step = 3
		-units = molality
	Results in	kmg_eut.jc%
			kmg_eut.jc&
			kmg_eut.jc@



	Sample 2: KMG_PER
	CAUTION: the results of the simulation are stored in files with
		-storage step = 3
		-units = molality
	Results in	kmg_per.jc%
			kmg_per.jc&
			kmg_per.jc@




HOW TO CREATE A SAMPLE FILE?
============================

Rather than introduce the simulation parameters by means of the
keyboard, users can create sample files that will be read by EQL.
(See the content of TEST_1.DAT or TEST_2.DAT)

Format of the sample file:

First record: 

	The names of the variables, separated by a comma (,) 
	or by one or more spaces (no spaces in the sample name).
	Recognized names are:

	label		sample name
	t		temperature in Deg. C
	ds		density (=1 if molalities)
	ph		pH
	alk		Alkalinity (in meq/l or meq/kg[H2O]) 
	na		total concentration for Na (in mmol/l or mmol/kg[H2O]) 
	k		K
	li		Li
	ca		Ca
	mg		Mg
	cl		Cl
	so4		SO4
	no3		NO3
	si		Si
	b		B

	Other variables (without spaces) can be present at any places in the data
 	file but they will be ignored by EQL.

	The order of variables is unimportant. The only rule is to set label as
	the first variable. 

	The following components must always be present (minimal set) :
	label, t, ds, ph, alk. 

	A component not analysed may be represented by: na  (or 0)
	A component analysed but not detected may be represented by: nd (or 0)
	In both cases the program converts na or nd to zero. No other alphanumeric
	expression is allowed. 

Analysis records:
	
	Each line/record contains the components values of one sample, separated 
	by a comma (,) or by one or more spaces. They must follow the same order 
	as defined in the first record. 

	Remark: If the concentrations are given in molalities, density has to be
	equaled to 1.



ERROR MESSAGES
==============

Most of them are self-explanatory. Here some informations about less
trivial errors.

In EQL:

"The equation set diverges: end of program"
When such a message is diplayed in the EQL initialisation program,
it is almost always due to inconsistent pH-alkalinity data or to 
a very bad electrical balance (no warning message is displayed if
the balance is too bad).     

"The system is in thermodynamic desequilibrium"
"The activity of water is constrained at different values"
"by more than one mineral assemblage"
Example: Let's consider the mineral assemblage:
	anhydrite (CaSO4)
	epsomite (MgSO4.7H2O)
	gypsum (CaSO4.2H2O)
	kieserite (MgSO4.H2O)
The activity of water is constrained by two assemblages:
	- anhydrite + gypsum
	- epsomite + kieserite
Each assemblage fixes a distinct water activity, which is impossible.

"System in thermodynamic desequilibrium: inconsistent mineral assemblage"
Let's consider the mineral assemblages:
Example 1: 
	anhydrite (CaSO4)
	glauberite (Na2Ca(SO4)2)
	thenardite (Na2SO4)
The solubility products of the three minerals would be constrained:
K(glauberite) = K(anhydrite) x K(thenardite), which is impossible.

Example 2: At t=25°C and log(PCO2) = -3.5:
	nahcolite (NaHCO3)
	natron (Na2CO3.10H2O)
cannot coexist. 	

In EVP:

"The equation system diverges: end of simulation"
The message occurs mostly in highly concentrated brines when small 
changes in concentrations produce high variations in activity 
coefficients. Reducing the convergence factor may sometimes fix the problem.
However, such cases are generally beyond the validity of Pitzer's equations.    



END CONDITIONS in EVP
=====================

If the program ends on an invariant point ("normal" end) the following 
messages are displayed: 
"invariant system / eutectic point" 
"invariant system / pseudo-eutectic point" 
"invariant system / peritectic point" 
"invariant system / pseudo-peritectic point" 

The term "pseudo" means that the solution has reached an invariant composition 
with a number of minerals in equilibrium lower than the maximum permitted by the 
phase rule.  

When the program ends "normally", the mineral file (sample_name.j[oc]%) is 
automatically compacted in order to suppress the minerals that have never
appeared during the simulation. When the simulation ends on an error message,
the user is asked whether he wants to compact the mineral file or not.
It is recommended to answer yes. If the program is stopped from the 
keyboard, no compaction is performed. 



HOW TO JOIN THE AUTHORS
=======================

François Risacher
I.R.D. - C.G.S.
1 rue Blessig
67084 Strasbourg Cedex
FRANCE
Tel: (33) (0)3 88 35 85 37
fax: (33) (0)3 88 36 72 35
email: risacher@illite.u-strasbg.fr

Alain CLEMENT
E.O.S.T. - C.G.S.
1 rue Blessig
67084 Strasbourg Cedex
FRANCE
Tel.: (33) (0)3 88 35 85 72
Fax:  (33) (0)3 88 36 72 35
email: aclement@illite.u-strasbg.fr
