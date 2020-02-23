import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import animation
import numpy as np
import tkinter as Tk
from tkinter import ttk
from functools import partial
from tkinter.colorchooser import askcolor
import math

G = 6.67e-11
class Body:
    def __init__(self, name, mass,xpos,ypos,xvel,yvel,index,bodylist,size=2,color="#000000"):
        self.mass = mass #define some properties for the body. For the simulation side of things, the most important ones are the mass, xpos, ypos, xvel and yvel.
        self.name = name
        self.size = size
        self.color = color
        self.xpos = xpos
        self.ypos = ypos
        self.xvel = xvel
        self.yvel = yvel
        self.index = index
        self.bodylist = bodylist
        self.derivative = [0,0,0,0]
        self.past_pos = [[],[]] #a 2d list that stores the body's past positions (to show off the shape of the orbit)
    
    def reset(self):
        global rel_var
        global rel_names
        
        new_x_pos = float(reset_x_pos.get())
        new_y_pos = float(reset_y_pos.get())
        new_x_vel = float(reset_x_vel.get())
        new_y_vel = float(reset_y_vel.get())
        rel_body = rel_names[str(rel_var.get())]
        if(rel_body!="COM" and rel_body!="Origin"):
            new_x_pos += rel_body.xpos
            new_y_pos += rel_body.ypos
            new_x_vel += rel_body.xvel
            new_y_vel += rel_body.yvel
        if(rel_body=="COM"):
            new_x_pos += com_pos[0]
            new_y_pos += com_pos[1]
            new_x_vel += com_vel[0]
            new_y_vel += com_vel[1]
        self.xpos = new_x_pos
        self.ypos = new_y_pos
        self.xvel = new_x_vel
        self.yvel = new_y_vel
        self.past_pos=[[],[]]
    
    def set_mass(self):
        self.mass = float(reset_mass.get())
        
    def make_circular_orbit(self): #use some quick maths to make the current body move on a circular orbit around the reference
        global rel_var
        global rel_names
        rel_body = rel_names[str(rel_var.get())] #get the reference body
        if(rel_body!="Origin" and rel_body !="COM"):
            rel_pos_x = self.xpos - rel_body.xpos 
            rel_pos_y = self.ypos - rel_body.ypos
            alt = np.sqrt(rel_pos_y**2 + rel_pos_x**2) #current dist from reference
        
            v_circular = np.sqrt(G*rel_body.mass/alt) #magnitude of the velocity of a circular orbit
        
            self.xvel = -v_circular*(rel_pos_y/alt)+rel_body.xvel #set the current body's velocity to that of a circular orbit RELATIVE to the reference
            self.yvel = v_circular*(rel_pos_x/alt)+rel_body.yvel
        '''elif(rel_body=="COM"):
            rel_pos_x = self.xpos - com_pos[0]
            rel_pos_y = self.ypos - com_pos[1]
            alt = np.sqrt(rel_pos_y**2 + rel_pos_x**2) #current dist from reference
        
            v_circular = np.sqrt(G*reduced_mass/alt) #magnitude of the velocity of a circular orbit
        
            self.xvel = -v_circular*(rel_pos_y/alt)+com_vel[0] #set the current body's velocity to that of a circular orbit RELATIVE to the reference
            self.yvel = v_circular*(rel_pos_x/alt)+com_vel[1]'''
        
        
    def dstate_dt(self): #create a function that returns the derivative of the body's position and velocity at any point
        self.radii = [0] * len(self.bodylist)
        self.x_accels = [0] *len(self.bodylist)
        self.y_accels = [0] *len(self.bodylist)
        for body in self.bodylist:
            if(body.index != self.index): #ensure that this body doesn't interact with itself
                self.radii[body.index] = np.sqrt((body.xpos-self.xpos)**2 + (body.ypos-self.ypos)**2) #get the distance between this body and all others
                
                self.x_accels[body.index] = -(G*body.mass*(self.xpos-body.xpos))/(self.radii[body.index]**3) #get the acceleration due to gravity from Newton's equation.
                
                self.y_accels[body.index] = -(G*body.mass*(self.ypos-body.ypos))/(self.radii[body.index]**3) #this is the derivative of the particle's velocity.
                
        
        state_deriv = [0,0,0,0]
        state_deriv[2] = sum(self.x_accels)
        state_deriv[3] = sum(self.y_accels)
        state_deriv[0] = self.xvel
        state_deriv[1] = self.yvel
        return state_deriv
    
    def step(self,d_t): #Uses the derivative given by dstate dt, along with the Euler Cromer method, to update the body's position and velocity.
        global rel_view_var
        global rel_names
        rel_view_body = rel_names[rel_view_var.get()]
        self.derivative = self.dstate_dt()
        self.xvel += d_t*(self.derivative[2])
        self.yvel += d_t*(self.derivative[3])
        self.xpos += d_t*self.xvel
        self.ypos += d_t*self.yvel
        if(rel_view_body!="Origin" and rel_view_body!="COM"):
            self.past_pos[0].append(self.xpos-rel_view_body.xpos) #this array stores past values for the object's position to allow them to be graphed
            self.past_pos[1].append(self.ypos-rel_view_body.ypos)
            
        elif(rel_view_body=="Origin"):
            self.past_pos[0].append(self.xpos) #this array stores past values for the object's position to allow them to be graphed
            self.past_pos[1].append(self.ypos)
        elif(rel_view_body=="COM"):
            self.past_pos[0].append(self.xpos-com_pos[0]) #this array stores past values for the object's position to allow them to be graphed
            self.past_pos[1].append(self.ypos-com_pos[1])
        if(len(self.past_pos[0])>15000): #delete sufficiently old values for the position to prevent messy orbit traces
                del(self.past_pos[0][0])
                del(self.past_pos[1][0])    
       


fig = plt.figure(figsize=[8,8])
ax = plt.axes()
ax.grid(which="both")

delta_t = 1 #timestep for simulation
time_accel = 30 #number of steps to skip for each frame of animation
time = 0
list_bodies = [] #list of all interacting bodies in simulation.

body1 = Body("Earth", 6e24,-2e8,0,0,-625,0,list_bodies,size=6,color="blue") #define some bodies with initial masses, positions and velocities to interact.
body2 = Body("Planet2", 5e24,2e8,0,0,750,1,list_bodies,size=5,color="green")
body3 = Body("Moon1", 3e21,0e8,0,-1000,400,2,list_bodies,size=3,color="gray")
plt.xlim(-5e8,5e8)
plt.ylim(-5e8,5e8)

list_bodies.append(body1) #add the bodies to the list
list_bodies.append(body2)
list_bodies.append(body3)



total_mass = sum(body.mass for body in list_bodies)


com_pos = [0,0] #define lists to store the position and velocity of the centre of mass
com_vel = [0,0]



lines = [] #variables that store the data to be plotted in each frame
dots = []

for body in list_bodies:
    lines.append(ax.plot([], [],ms=0.5,color=body.color)[0])
    dots.append(ax.plot([],[],"o",ms=body.size,color=body.color)[0])




paused = False
def pause(): # pause or unpause the simulation
    global paused
    paused = not paused
    
def current_reset(): #reset the current body, giving it new positions and velocities
    global body_names
    body_names[name_var.get()].reset()

def current_set_mass(): #set the mass of the current body
    global body_names
    body_names[name_var.get()].set_mass()
    
def fill_current_position(clear): #if clear is true, clear the fields. Otherwise set them to the object's current position
    global body_names
    global rel_var
    global rel_names
    rel_body = rel_names[rel_var.get()]
    curr_x_pos = round(body_names[name_var.get()].xpos,2)
    curr_y_pos = round(body_names[name_var.get()].ypos,2)
    if(rel_body!="Origin" and rel_body!="COM"):
        curr_x_pos = round(body_names[name_var.get()].xpos-rel_body.xpos,2)
        curr_y_pos = round(body_names[name_var.get()].ypos-rel_body.ypos,2)
    if(rel_body=="COM"):
        curr_x_pos = round(body_names[name_var.get()].xpos-com_pos[0],2)
        curr_y_pos = round(body_names[name_var.get()].ypos-com_pos[1],2)
        
    reset_x_pos.delete(0,Tk.END)
    reset_y_pos.delete(0,Tk.END)
    if(not clear):
        reset_x_pos.insert(0,curr_x_pos)
        reset_y_pos.insert(0,curr_y_pos)
    
def fill_current_velocity(clear): #same as the function to fill positions. Fills relative to reference body
    global body_names
    global rel_var
    global rel_names
    rel_body = rel_names[rel_var.get()]
    
    curr_x_vel = round(body_names[name_var.get()].xvel,2)
    curr_y_vel = round(body_names[name_var.get()].yvel,2)
    
    if(rel_body!="Origin" and rel_body!="COM"):
        curr_x_vel = round(body_names[name_var.get()].xvel-rel_body.xvel,2)
        curr_y_vel = round(body_names[name_var.get()].yvel-rel_body.xvel,2)
    if(rel_body=="COM"):
        curr_x_vel = round(body_names[name_var.get()].xvel-com_vel[0],2)
        curr_y_vel = round(body_names[name_var.get()].yvel-com_vel[1],2)
    
    reset_x_vel.delete(0,Tk.END)
    reset_y_vel.delete(0,Tk.END)
    if(not clear):
        reset_x_vel.insert(0,curr_x_vel)
        reset_y_vel.insert(0,curr_y_vel)
        
def set_color(): #change the body's colour using tkinter.colorchooser
    global body_names
    
    curr_body =  body_names[name_var.get()]
    new_color = askcolor()
    if(new_color[1] == None):
        new_color[1]= curr_body.color
    curr_body.color = new_color[1]
    lines[curr_body.index].set_color(new_color[1])
    dots[curr_body.index].set_color(new_color[1])
    
def set_circular_orbit(): #call the class method which moves the current body onto a circular orbit relative to the reference body (the altitude is the same)
    global body_names
    
    curr_body =  body_names[name_var.get()]
    curr_body.make_circular_orbit()

def delete_current_body(name):
    global body_names
    global list_bodies
    
def clear_trails(*args,widget=None):
    for body in list_bodies:
        body.past_pos=[[],[]]
    
def timestep(speed_up,accelerate): #function to change timestep when buttons are pressed
    #first argument is bool, tells whether to increase the variable or decrease it
    #second argument is bool, tells which variable to change - the time step used to simulate, or the acceleration (no. of steps taken per frame of animation)
    global delta_t
    global time_accel
    change_variable = delta_t
    if(accelerate):
        change_variable = time_accel
    if(speed_up):
        if(change_variable<1):
            change_variable+=0.1
            change_variable = round(change_variable,1)
        elif(change_variable>=1 and change_variable<10):
            change_variable += 1
        elif(change_variable>=10 and change_variable<100):
            change_variable += 5
        elif(change_variable>=100):
            change_variable += 20
    else:
        if(change_variable<=1 and change_variable>0.1):
            change_variable-=0.1
            change_variable = round(change_variable,1)
        elif(change_variable<=10 and change_variable>1):
            change_variable -= 1
        elif(change_variable>10 and change_variable <= 100):
            change_variable -= 5
        elif(change_variable>100):
            change_variable -= 20
    if(accelerate):
        time_accel = round(change_variable)
    else:
        delta_t = change_variable

    

root = Tk.Tk()
root.geometry("1100x1000") #set window size
root.title("Principia")
frame = Tk.Frame(root)

speed_var = Tk.StringVar(root)
speed_label = ttk.Label(root,textvariable=speed_var) #label to display speed of current object

mass_var = Tk.StringVar(root)
mass_label = ttk.Label(root,textvariable=mass_var) #label to display current mass

time_label = ttk.Label(root,text="Time: ") #label displaying time

body_drop_label = ttk.Label(root,text="Current Body") #dropdown menus for the current body being observed, the reference body, and the current VIEW reference.
rel_drop_label = ttk.Label(root,text="Reference")
rel_view_label = ttk.Label(root,text="Relative View")

alt_var = Tk.StringVar(root)
alt_label = ttk.Label(root,textvariable=alt_var) #label to display altitude of current object relative to reference.


canvas = FigureCanvasTkAgg(fig, master=root)
toolbar = NavigationToolbar2Tk( canvas, frame )

name_var = Tk.StringVar(root)
rel_var = Tk.StringVar(root)

rel_view_var = Tk.StringVar(root)
rel_view_var.trace("w",partial(clear_trails,widget=rel_view_var))

reset_var = Tk.StringVar(root)

body_names ={body.name : body for body in list_bodies} #dictionary with keys of body names and values of the corresponding body
rel_names = body_names.copy()
rel_names.update({"Origin":"Origin"}) #dictionary for relative names to appear in relative dropdown
rel_names.update({"Centre of Mass":"COM"}) #dictionary for relative names to appear in relative dropdown

body_menu = ttk.OptionMenu(root,name_var,body1.name,*body_names) #dropdowns for the current body, reference body, and view reference.
rel_menu  = ttk.OptionMenu(root,rel_var,"Origin",*rel_names)
rel_view_menu  = ttk.OptionMenu(root,rel_view_var,"Origin",*rel_names)


new_x_pos = 0
new_y_pos = 0
new_x_vel = 0
new_y_vel = 0


reset_x_pos = ttk.Entry(root) #Entries to reset the body's position, velocity and mass.
reset_y_pos = ttk.Entry(root)
reset_x_vel = ttk.Entry(root)
reset_y_vel = ttk.Entry(root)

reset_mass = ttk.Entry(root)

reset_pos_label = ttk.Label(root,text="New X/Y Positions (m)")
reset_vel_label = ttk.Label(root,text="New X/Y Velocities (m/s)")


reset_x_pos.insert(0,"1e7")
reset_y_pos.insert(0,"1e7")
reset_x_vel.insert(0,"1")
reset_y_vel.insert(0,"1")

pause_button = ttk.Button(root,text="Pause",command = partial(pause)) 


timestep_label = ttk.Label(root,text="Timestep: ")
timestep_plus = ttk.Button(root,text="+",command = partial(timestep,True,False))
timestep_minus = ttk.Button(root,text="-",command = partial(timestep,False,False))

time_accel_label = ttk.Label(root,text="Time Acceleration: ")
time_accel_plus = ttk.Button(root,text="+",command = partial(timestep,True,True))
time_accel_minus = ttk.Button(root,text="-",command = partial(timestep,False,True))


set_current_position = ttk.Button(root,text="Current",command=partial(fill_current_position,False))
set_current_velocity = ttk.Button(root,text="Current",command=partial(fill_current_velocity,False))
clear_positions = ttk.Button(root,text="Clear",command=partial(fill_current_position,True))
clear_velocities = ttk.Button(root,text="Clear",command=partial(fill_current_velocity,True))
reset_button = ttk.Button(root,text="Set Current Body State",command=partial(current_reset))

reset_mass_button = ttk.Button(root,text="Change Mass",command=partial(current_set_mass))

color_button = ttk.Button(root,text = "Change Colour", command = partial(set_color))

circular_orbit_button = ttk.Button(root,text="Make Circular Orbit around Reference", command = partial(set_circular_orbit))

#place all the widgets listed above.

rel_drop_label.place(x=125,y=80) 
body_drop_label.place(x=10,y=80)

body_menu.place(x=10,y=100,height=25,width=110)
rel_menu.place(x=125,y=100,height=25,width=110)

pause_button.place(x=10,y=150)

alt_label.place(x=10,y=220)
speed_label.place(x=10,y=200)
mass_label.place(x=10,y=180)
time_label.place(x=10,y=240)

rel_view_label.place(x=140,y=150)
rel_view_menu.place(x=140,y=170)

timestep_plus.place(x=50,y=300,height = 30, width=30)
timestep_minus.place(x=10,y=300,height=30, width = 30)
timestep_label.place(x=10,y=280)

time_accel_plus.place(x=50,y=360,height = 30, width=30)
time_accel_minus.place(x=10,y=360,height=30, width = 30)
time_accel_label.place(x=10,y=340)

reset_pos_label.place(x=10,y=400)
reset_x_pos.place(x=10,y=420,height = 30, width = 60)
reset_y_pos.place(x=70,y=420,height = 30, width = 60)
set_current_position.place(x=140,y=415)
clear_positions.place(x=140,y=437)

reset_vel_label.place(x=10,y=470)
reset_x_vel.place(x=10,y=490,height = 30, width = 60)
reset_y_vel.place(x=70,y=490,height = 30, width = 60)
set_current_velocity.place(x=140,y=485)
clear_velocities.place(x=140,y=507)

reset_button.place(x=10,y=530,height = 30, width = 150)

circular_orbit_button.place(x=10,y=570,height=30,width=230)

reset_mass.place(x=10,y=650,height=20,width=80)
reset_mass_button.place(x=100,y=650)

color_button.place(x=10,y=700)


canvas.get_tk_widget().place(x=250,y=100)
frame.place(x=250,y=70)


def init(): #called when the animation begins.

    ax.set_aspect('equal')

    return tuple(lines)+tuple(dots)

def animate(n):
    global time
    if(not paused): #if the simulation is running, update the bodies
        for body in list_bodies:
            for i in range(0,time_accel):
                body.step(delta_t)
        time += (delta_t*time_accel)
    
    curr_body = body_names[str(name_var.get())] #get the current body (given by the "current" dropdown) by looking up the dictionary
    rel_view_body = rel_names[rel_view_var.get()] #get the body currently being used as a view reference.
    #Note that rel_view_body is not the same as rel_body; rel_body is the body which the current body's position and velocity is displayed relative to.
    #rel_view_body is the body that the simulation is VIEWED with respect to.
    
    for body in list_bodies: #set the line and dot data, whether the simulation is running or not. This way the bodies reset properly even when the sim is paused
        lines[body.index].set_data(body.past_pos[0],body.past_pos[1])
        if(rel_view_body!="Origin" and rel_view_body!="COM"):
            dots[body.index].set_data(body.xpos-rel_view_body.xpos,body.ypos-rel_view_body.ypos)
        elif(rel_view_body=="Origin"):
            dots[body.index].set_data(body.xpos,body.ypos)
        elif(rel_view_body=="COM"):
            dots[body.index].set_data(body.xpos-com_pos[0],body.ypos-com_pos[1])
    com_pos[0] = (sum(body.mass*body.xpos for body in list_bodies))/total_mass #get the state of the center of mass
    com_pos[1] = (sum(body.mass*body.ypos for body in list_bodies))/total_mass
    com_vel[0] = (sum(body.mass*body.xvel for body in list_bodies))/total_mass
    com_vel[1] = (sum(body.mass*body.yvel for body in list_bodies))/total_mass
    
    mass = curr_body.mass
    
    if(rel_var.get() == "Origin"): #if the reference frame is wrt origin, just square the body's velocity and position components
        speed = np.sqrt(curr_body.xvel**2 + curr_body.yvel**2)
        alt = np.sqrt(curr_body.xpos**2 + curr_body.ypos**2)
        
    elif(rel_var.get() == "Centre of Mass"): #ditto for the center of mass.
        speed = np.sqrt((curr_body.xvel-com_vel[0])**2 + (curr_body.yvel-com_vel[1])**2)
        alt = np.sqrt((curr_body.xpos-com_pos[0])**2 +(curr_body.ypos-com_pos[1])**2)
    else:
        rel_body = rel_names[str(rel_var.get())] #otherwise, output the speed and altitude relative to the reference body.
        speed = np.sqrt((curr_body.xvel-rel_body.xvel)**2 +(curr_body.yvel-rel_body.yvel)**2)
        alt = np.sqrt((curr_body.xpos-rel_body.xpos)**2 + (curr_body.ypos-rel_body.ypos)**2)
        
    speed_var.set("Speed: " + str(round(speed,2))+" m/s") #update speed, altitude, timestep, time_accel, and time for the GUI readouts
    alt_var.set("Altitude: " +str(round((alt/1000),2))+" km")
    mass_var.set("Mass: "+str(round(mass,3))+" kg")
    
    timestep_label.configure(text="Timestep: "+str(delta_t)+" s")
    time_accel_label.configure(text="Time Acceleration: "+str(time_accel*100)+"x")
    time_label.configure(text="Time: "+ str(round(time/3600))+" hrs")
    
    return tuple(lines)+tuple(dots) #return the data to be graphed for the animation
    
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2000, interval=10, blit=True, repeat=True) #animate the matplotlib plot

canvas.draw() #update the tkinter canvas with the matplotlib animation

Tk.mainloop()
