I will write some compiling notes and documentation here

First of all make your own branch of the code so we can work on them independently.

Compiling:

For your installation you probably want to do:

>cp Makefile_AndrewLinux Makefile_TylerLinux

and then edit it so that NAGCDIR looks like what mine looks like, but for your distribution of NAG

The intel64/libi libraries are included as per a guide i found on the web, though I don't know what they are. Theyre specific to the linux version of NAG. Once you've editted it then

>cp Makefile_TylerLinux Makefile
>make

should give you the executable md_out. 

md_test.c contains the main() which I'm working with currently. You can create your own working function in that file and call it from main. Look at simple_test() and see how that creates a Polymer, simulation parameters, and feeds it to simulation(). simulation() then calls update() once every time step, which does everything in the following order:

-compute forces on all particles (6πµv for those particles whose motion v is perscribed)
-calculate a noise vector (zero for those particles whose motion is perscribed)
-calcuate D tensor
-calculate D*f
-decompose D into B
-calculate B*dw (noise)
-update particle positions

That is basically most there is to it.

You should probably aim to compute the angular potential in computeForces(). You will spot the space I've left for your 2D RP tensor. BC's Im not too sure about. If they involve a lot of structural change then let me know because it will make it harder to merge branches later. Hope it goes well! Let me know if you have any problems
