'''
create an empty particle stack
read open particle set and get at list of new stacks
if new stacks, iterate through the list of new stacks and insertfunctionstep to move each stack
	get the following particle stack, if there is none, go back to sleep
	check if there is enough disk space for the destination stacks
		if not enough space, stop here
	put the particles in the stack into a loading port
	spawn a new copy operation for the next stack
	when copy finished, make the new particle stack available
else
	sleep
go back to 2
'''
