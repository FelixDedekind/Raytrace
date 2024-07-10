import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy.interpolate import interp1d,UnivariateSpline
from scipy.optimize import root_scalar
from scipy.misc import derivative

#constants:

n_water = 1.33
n_air = 1
n_glass = 1.4
# Sample data points for the curve



class beam:
    def __init__(self, x0, y0,direction,colour):
        self.path_x = [x0]
        self.path_y = [y0]
        self.direction = direction
        self.colour = colour
        self.d = 0

    def line(self,x0):
        x = self.path_x[-1]
        y = self.path_y[-1]
        d = y - self.direction * x
        self.d = d
        return self.direction*x0 + d


class lens:
    def __init__(self,x_vals,y_vals):
        self.x_vals = x_vals
        self.y_vals = y_vals
        self.shape_curve =interp1d(x_vals,y_vals, kind='cubic', fill_value="extrapolate") 
        self.spline = UnivariateSpline(x_vals,y_vals,s=0)
        self.derivative = self.spline.derivative()





class setup:
    def __init__(self,lens,beams,dl,do,radius,num_pixels,glass_width, colour_func,target_size):
        self.lens = lens
        self.beams = beams
        self.dl = dl
        self.do = do
        self.lens_radius= radius
        self.num_pixels = num_pixels
        self.glass_width = glass_width
        self.colour_func = colour_func
        self.target_size = target_size
    
    def generate_viewer_beams(self, n_beams):
        res = []
        for i in range(n_beams):
            new_beam = beam(0,0,self.lens_radius/self.dl/n_beams *(i),-1)
            res.append(new_beam)
        self.beams = res

    def generate_point_beams(self,dy, n_beams):
        res = []
        for i in range(n_beams):
            new_beam = beam(0,dy,0.1*self.lens_radius/self.dl/n_beams * 2*(i+1-n_beams/2),-1)
            res.append(new_beam)
        self.beams = res
    
    
    def generate_parallel_beams(self,n_beams):
        res = []
        for i in range(n_beams):
            new_beam = beam(0,radius/n_beams*i,0,-1)
            res.append(new_beam)
        self.beams = res
        
    


        


def plot_setup(setup,plot_focus,n):
    fig,ax = plt.subplots()
    ax.set_xlim((0,setup.dl+setup.do*1.1))

    #plot lens:
    xlens = setup.lens.x_vals
    ylens = setup.lens.y_vals
    xlens = np.append(xlens,np.flip(xlens))
    ylens = np.append(ylens,-np.flip(ylens))
    ax.plot(xlens+setup.dl,ylens,label='lens',c='black')
    #ax.plot(xlens+setup.dl,setup.lens.shape_curve(xlens))

    #plot glass:
    xs = np.array([0,0,setup.glass_width,setup.glass_width,0])+setup.dl+np.max(setup.lens.x_vals)
    ys = np.array([1,-1,-1,1,1])*setup.lens_radius
    ax.plot(xs,ys,c='black')
    


    #plot beams:
    d_sum = 0
    k_sum = 0
    beam_i = 0

    for beam in setup.beams:
        if beam_i % n != 0: 
            beam_i+=1
            continue
        beam_i +=1
        d_sum += beam.d
        k_sum += beam.direction
        colour = 'black'
        if beam.colour == 1:
            colour = 'red'
        ax.plot(beam.path_x,beam.path_y,c=colour)
        #ax.scatter(beam.path_x,beam.path_y)
        ax.plot(beam.path_x,-np.array(beam.path_y), c=colour)
        #ax.scatter(beam.path_x,-np.array(beam.path_y))
        #xs = np.linspace(beam.path_x[-1],setup.dl+setup.do)
        #ax.plot(xs,beam.line(xs))
        #ax.plot(xs,-beam.line(xs))

    
    if plot_focus: 
        x = -d_sum/k_sum
        xs = [x,x]
        ys = [-setup.lens_radius,setup.lens_radius]
        ax.plot(xs,ys)
        print(f"Die Brennweite ist {x-setup.dl-np.max(xlens)}")
    
    #plot target
    xs = np.array([1,1])*(setup.dl+setup.do)
    ys = np.array([-1,1])
    plt.plot(xs,ys*setup.target_size,color = 'red')
    plt.plot(xs,ys*setup.lens_radius*5,color = 'black')
    
    
    

    plt.show()

def read_lens_shape(radius,path):
    df = pd.read_csv(path)
    xspace = df['x'].to_numpy()
    yspace = df['y'].to_numpy()


    xspace = xspace * -1
    #yspace = yspace * -1

    xspace -= np.min(xspace)
    yspace -= np.min(yspace)

    scale = radius/np.max(np.abs(xspace))

    xspace *= scale
    yspace *= scale

    sorted_indices_y = np.argsort(xspace)
    sorted_indices_x = np.argsort(yspace)


    xs = yspace[sorted_indices_x]
    ys = xspace[sorted_indices_y]
    
    #unique_xs, unique_is = np.unique(xs, return_index=True)
    #unique_ys = ys[unique_is]

    
    #plt.plot(xs,ys)
    #plt.show()



    return lens(xs,ys)

def trace_beams(setup):
    for beam in setup.beams:
        #find point where lens is hit
        func = lambda x: -beam.line(x) + setup.lens.shape_curve(x-setup.dl)
        min = np.min(setup.lens.x_vals)
        max = np.max(setup.lens.x_vals)
        domain = (setup.dl,max+setup.dl)
        try:
            result = root_scalar(func,bracket = domain)
            x_intersection = result.root
            #print(f"x intersection is {x_intersection} the y value is {setup.lens.shape_curve(x_intersection-setup.dl)} remainder is {func(x_intersection)}")
            beam.path_x.append(x_intersection)
            beam.path_y.append(setup.lens.shape_curve(x_intersection-setup.dl))
        except Exception as e:
            xs = np.linspace(domain[0],domain[1],1000)
            ys = func(xs)
            #plt.plot(xs,ys)
            #plt.plot(xs,setup.lens.shape_curve(xs-setup.dl))
            #plt.show()
            print(f"An error occured: {e}")
            continue

        #refract when hitting the lens 
        #derivative_at_beam = derivative(setup.lens.shape_curve,beam.path_x[-1]-setup.dl,dx = 1e-12)
        derivative_at_beam = setup.lens.derivative(beam.path_x[-1]-setup.dl)

        angle_beam = np.arctan(beam.direction)
        angle_lens = np.arctan(derivative_at_beam)

        alpha = np.pi/2 - angle_lens + angle_beam

        beta = np.arcsin(n_air/n_water*np.sin(alpha))
        gamma = angle_lens-np.pi/2+beta

        beam.direction = np.tan(gamma)

        #print(f"y {setup.lens.shape_curve(x_intersection-setup.dl)} derivative: {derivative_at_beam} angle lens: {angle_lens} angle beam: {angle_beam} alpha: {alpha} beta: {beta}")
        xs = np.linspace(beam.path_x[-1],setup.dl+setup.do,1000)
        ys = beam.line(xs)
        #plt.plot(xs,ys)
        #plt.show()

        #refract at water glass
        xnew = setup.dl+np.max(setup.lens.x_vals)
        ynew = beam.line(xnew)

        beam.path_x.append(xnew)
        beam.path_y.append(ynew)

        alpha = np.arctan(beam.direction)
        beta = np.arcsin(n_water/n_glass*np.sin(alpha))

        beam.direction = np.tan(beta)
        #refract at glass air
        xnew = setup.glass_width+beam.path_x[-1]
        ynew = beam.line(xnew)

        beam.path_x.append(xnew)
        beam.path_y.append(ynew)

        alpha =np.arctan(beam.direction)
        beta = np.arcsin(n_glass/n_air*np.sin(alpha))

        beam.direction = np.tan(beta)
        #add screen point

        xnew = setup.dl + setup.do
        ynew = beam.line(xnew)

        beam.path_y.append(ynew)
        beam.path_x.append(xnew)
        #print(f"end of x {beam.path_x[-1]} and y {beam.path_y[-1]} and calculated {beam.line(beam.path_x[-1])} and direction {beam.direction}")


def colour_beams(setup):
    for beam in setup.beams:
        beam.colour = setup.colour_func(beam.path_y[-1],setup.target_size)
'''
def colour_pixels(setup):
    pixels = np.zeros(setup.num_pixels)
    n_rays = np.zeros(setup.num_pixels)
    for beam in setup.beams:
        if len(beam.path_x) < 2:
            print("wtf")
            continue
        index = beam.path_y[1]/setup.lens_radius*setup.num_pixels
        print(f"index = {index} y value is {beam.path_y[1]} the radius is {setup.lens_radius} colour is: {beam.colour}")
        index = int(index)

        n_rays[index]=n_rays[index]+1.
        pixels[index]=pixels[index]+beam.colour
    pixels_normalized = pixels
    for i in range(setup.num_pixels):
        if(n_rays[i]!=0):
            pixels_normalized[i] /= n_rays[i]


    #print(pixels)
    return n_rays,pixels,pixels_normalized
'''
def colour_pixels(setup):
    # Initialize arrays for pixel colors and ray counts
    pixels = np.zeros(setup.num_pixels)
    n_rays = np.zeros(setup.num_pixels)
    
    for beam in setup.beams:
        # Check if the beam path has at least two points
        if len(beam.path_x) < 2:
            print("Beam path is too short")
            continue
        
        # Compute the index for the beam's interaction with the lens
        index = beam.path_y[1] / setup.lens_radius * setup.num_pixels
        #print(f"Index = {index:.2f}, y value = {beam.path_y[1]:.2f}, radius = {setup.lens_radius:.2f}, colour = {beam.colour}")
        
        # Convert to integer index
        index = int(index)
        
        # Ensure index is within valid range
        if 0 <= index < setup.num_pixels:
            n_rays[index] += 1
            pixels[index] += beam.colour
        else:
            print(f"Index {index} is out of bounds")
    
    # Normalize pixel colors
    pixels_normalized = np.zeros_like(pixels)
    nonzero_indices = n_rays != 0
    pixels_normalized[nonzero_indices] = pixels[nonzero_indices] / n_rays[nonzero_indices]
    
    return n_rays, pixels, pixels_normalized

        

mylens = read_lens_shape(5e-3,"shape/out/200-smoothed.csv")

radius = 0.005
xs = np.linspace(0,radius,1000)
ys = np.sqrt(radius**2-xs**2)
xs = radius-xs
#plt.plot(xs,ys)
#jplt.xlim((0,2*radius))
#plt.ylim((0,2*radius))
#plt.show()
#mycircularlens = lens(xs,ys)


ys = 10*xs**2
#plt.plot(ys,xs)
#plt.show()
#myparaboliclens = lens(ys,xs)




def colour_func(y,target_size):
    if np.abs(y) < target_size:
        return 1.
    return 0

mysetup = setup(
    lens = mylens, 
    beams = 0,
    dl=0.1,
    do=0.015,
    radius=0.005,
    num_pixels=100,
    glass_width=1e-3,
    colour_func = colour_func,
    target_size = 1e-3
    )
n_beams = 10000
mysetup.generate_viewer_beams(n_beams)
#print(len(mysetup.beams))
#mysetup.generate_parallel_beams(50)
trace_beams(mysetup)
colour_beams(mysetup)



n_rays,pixels,pixels_normalized = colour_pixels(mysetup)

xs = np.linspace(0,mysetup.lens_radius,mysetup.num_pixels)

data = {
    'xs': xs,
    'n_rays': n_rays,
    'pixels_normalized': pixels_normalized
}

df = pd.DataFrame(data)

df.to_csv("out.csv")

plt.plot(xs,pixels)
plt.scatter(xs,pixels)
plt.show()

first_drop_index = 0
threshold = 0.9

while(pixels[first_drop_index]>threshold):
    first_drop_index += 1

cutoff_distance = xs[first_drop_index]

print(f"cutoff_distance = {cutoff_distance}")





plot_setup(mysetup,True,1000)
