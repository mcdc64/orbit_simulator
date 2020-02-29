Cian McDonnell 2020
----- Orbit Simulator 1.0-----

A simulator that calculates the motions of objects which affect each other gravitationally.
Run the .py file to start the simulation. Requires tkinter, matplotlib, and numpy.

It is not currently possible to add new bodies while the simulation is running.
Instead, a new scenario can be created by changing the bodies which are hardcoded.
To change the scenario, edit the "Body" declarations on line 136 to whatever initial conditions you like.
Don't forget to change the axis limits, and to append whatever bodies you have to "list_bodies". (Lines 139-144).

To observe one of the scenarios included, copy the code you want from scenarios.txt and replace lines 136-144 in orbit_sim_gui.py
with that code.
Then run the simulation.
