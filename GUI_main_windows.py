import tkinter as tk
from tkinter import ttk, messagebox
import subprocess
import os
import xarray as xr
import numpy as np
import xgrads
import matplotlib.pyplot as plt
from PIL import Image, ImageTk

# Global variable for result_label
result_label = None

class FortranCompilerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Simple Cold Pool Hires GUI Interactive Model - Compiler")

        # Set the background image
        script_directory = os.path.dirname(__file__)
        background_image_path = os.path.join(script_directory, "background.png")

        try:
            # Load the background image using Pillow
            image = Image.open(background_image_path)
            self.background_image = ImageTk.PhotoImage(image)

            # Create a Canvas to display the background image
            canvas = tk.Canvas(self.root, width=image.width, height=image.height)
            canvas.pack()

            # Place the background image on the canvas
            canvas.create_image(0, 0, anchor=tk.NW, image=self.background_image)

            # Create a transparent frame for the title and button
            title_frame = ttk.Frame(self.root, padding=(10, 5))
            title_frame.place(relx=0.5, rely=0.5, anchor="center")

            # Add a label for the subtitle
            subtitle_label = ttk.Label(title_frame, text="Fortran Code Compilation", font=("Helvetica", 14))
            subtitle_label.grid(row=0, column=0, padx=10, pady=5)

            # Button to compile Fortran code
            self.compile_button = tk.Button(title_frame, text="Compile Fortran Code", width=30, height=15, bg="blue", fg="black", command=self.compile_fortran)
            self.compile_button.grid(row=2, column=0, pady=20)

        except Exception as e:
            messagebox.showerror("Error", f"Error loading background image: {e}")

      #  self.create_widgets()

    def create_widgets(self):
        # Button to compile Fortran code
        self.compile_button = tk.Button(self.root, text="Compile Fortran Code", width=30, height=15, bg="blue", fg="black", command=self.compile_fortran)
        self.compile_button.grid(row=2, column=0, pady=20)

    def compile_fortran(self):
        self.compile_button["state"] = "disable"
        
        makefile_config = """
FC = gfortran -O2 -fconvert=big-endian
FCFLAGS = 
"""

        makefile_config += f"""
SOURCES=src/DTDM_model.f src/blktri.f

OBJS= $(SOURCES:.f=.o)

dtdm: clean ; 
\t$(FC) $(FCFLAGS) -o $@ $(SOURCES) 

clean:
\t/bin/rm -f dtdm src/*.o
"""

        # Writing makefile configuration to a temporary file
        with open("temp_makefile", "w") as temp_makefile:
            temp_makefile.write(makefile_config)

        # Running the makefile
        try:
            subprocess.run(["make", "-f", "temp_makefile"])
            messagebox.showinfo("Compilation Complete", "Fortran code compiled successfully.")
        except Exception as e:
            messagebox.showerror("Compilation Error", f"Error compiling Fortran code: {e}")
        finally:
            # Removing temporary makefile
            os.remove("temp_makefile")

        # Enable compile button
        self.compile_button["state"] = "normal"

        # Close the main window
        self.root.destroy()

        # Create and display the second GUI window
        show_second_gui()

def show_second_gui():
    global result_label  # Access the global variable
    # Create the second GUI window
    app = tk.Tk()
    app.title("Simple Cold Pool Hires GUI Interactive Model - User Data Entry")
    
    app.geometry("1366x768")
  
    # Add input fields and buttons to the transparent frame
    title_frame = ttk.Frame(app, padding=(10, 5))
    title_frame.place(relx=0.5, rely=0.5, anchor="center")
    
    # Add input fields
    nx_label = ttk.Label(app, text="nx = ", font=("Helvetica", 16))
    nx_label.grid(row=0, column=0, padx=80, pady=40, sticky="e")
    nx_entry = ttk.Entry(app)
    nx_entry.grid(row=0, column=1, padx=80, pady=40, sticky ="e")

    nz_label = ttk.Label(app, text="nz = ",font=("Helvetica", 16))
    nz_label.grid(row=1, column=0, padx=80, pady=40, sticky="e")
    nz_entry = ttk.Entry(app)
    nz_entry.grid(row=1, column=1, padx=80, pady=40)

    dx_label = ttk.Label(app, text="dx =",font=("Helvetica", 16))
    dx_label.grid(row=2, column=0, padx=80, pady=40, sticky="e")
    dx_entry = ttk.Entry(app)
    dx_entry.grid(row=2, column=1, padx=80, pady=40)

    dz_label = ttk.Label(app, text="dz = ",font=("Helvetica", 16))
    dz_label.grid(row=3, column=0, padx=80, pady=40, sticky="e")
    dz_entry = ttk.Entry(app)
    dz_entry.grid(row=3, column=1, padx=80, pady=40)

    dt_label = ttk.Label(app, text="dt = ",font=("Helvetica", 16))
    dt_label.grid(row=4, column=0, padx=80, pady=40, sticky="e")
    dt_entry = ttk.Entry(app)
    dt_entry.grid(row=4, column=1, padx=80, pady=40)

    tim_label = ttk.Label(app, text="timend =",font=("Helvetica", 16))
    tim_label.grid(row=0, column=2, padx=100, pady=40, sticky="e")
    tim_entry = ttk.Entry(app)
    tim_entry.grid(row=0, column=3, padx=100, pady=40)

    plot_label = ttk.Label(app, text="plot =",font=("Helvetica", 16))
    plot_label.grid(row=1, column=2, padx=100, pady=40, sticky="e")
    plot_entry = ttk.Entry(app)
    plot_entry.grid(row=1, column=3, padx=100, pady=40)

    icz_label = ttk.Label(app, text="ICOOLZONE (input 0/1/2) = ",font=("Helvetica", 16))
    icz_label.grid(row=2, column=2, padx=100, pady=40, sticky="e")
    icz_entry = ttk.Entry(app)
    icz_entry.grid(row=2, column=3, padx=100, pady=40)

    ips_label = ttk.Label(app, text="IPRESSURE (input 0/1) = ",font=("Helvetica", 16))
    ips_label.grid(row=3, column=2, padx=100, pady=40, sticky="e")
    ips_entry = ttk.Entry(app)
    ips_entry.grid(row=3, column=3, padx=100, pady=40)

    ims_label = ttk.Label(app, text="IMOIST (input 0/1) =",font=("Helvetica", 16))
    ims_label.grid(row=4, column=2, padx=100, pady=40, sticky="e")
    ims_entry = ttk.Entry(app)
    ims_entry.grid(row=4, column=3, padx=100, pady=40)

    # Add a button to save to file
    save_button = ttk.Button(app, text="Save to File", command=lambda:save_to_file(nx_entry,nz_entry,dx_entry,dz_entry,dt_entry,tim_entry,plot_entry,icz_entry,ips_entry,ims_entry))
    save_button.grid(row=5, column=1, columnspan=2, pady=15)

    # Add a button to run the model
    run_model_button = ttk.Button(app, text="Run Model", command=lambda:run_model(True))
    run_model_button.grid(row=6, column=1, columnspan=4, pady=15)

    # Add a button to view info
    info_button = ttk.Button(app, text="Show Info", command=info)
    info_button.grid(row=5, column=2, columnspan=2, pady=15)

    # Display the result label
    result_label = ttk.Label(app, text="",font=("Helvetica", 14))
    result_label.grid(row=8, column=2, columnspan=2, pady=15)

    quit_button = ttk.Button(app, text="Quit", command=app.quit)
    quit_button.grid(row=7, column=1, columnspan=4, pady=15)
    
    

    # Start the GUI application
    app.mainloop()

def save_to_file(nx_entry,nz_entry,dx_entry,dz_entry,dt_entry,tim_entry,plot_entry,icz_entry,ips_entry,ims_entry):
    global result_label  # Access the global variable
    #user input:
    nx = nx_entry.get()
    nz = nz_entry.get()
    dx = dx_entry.get()
    dz = dz_entry.get()
    dt = dt_entry.get()
    tim = tim_entry.get()
    plot = plot_entry.get()
    icz = icz_entry.get()
    ips = ips_entry.get()
    ims = ims_entry.get()

    #create string
    data = f"""
 &experiment
    casename = 'coldpool.hires',
 /
 &grid_run
    nx = {nx},
    nz = {nz},
    dx = {dz}.,
    dz = {dz}.,
    dt = {dt},
    timend = {tim}.,
    plot = {plot}.,
 /
 &framework
    ipressure = {ips},
    ianelastic = 0,
    csnd = 50.,
    imoist = {ims},
    fcoef = 0.,
 /
 &numerics
    byteswap = 1,
    cstar = 30.,
    dkx = 50.,
    dkz = 10.,
    eps = 0.005,
 /
 &environ
    bvpbl = 0.000,
    pbld = 4000.,
    bvtropo = 0.00,
    tropo = 12000.,
    bvstrat = 0.00,
    psurf = 965.,
    usurf = 0.,
    shear1 = 0.,
    sdepth1 = 3000.,
    shear2 = 0.,
    sdepth2 = 1500.,
    shear3 = 0.,
    rhsurf = 80,
    rhvalue1 = 90,
    rhheight1 = 1000.,
    rhvalue2 = 85,
    rhheight2 = 4000.,
    rhvalue3 = 25,
    rhheight3 = 8000.,
 /
 &thermal
    ithermal = 0,
    delt = 3.,
    radx = 8000.,
    radz = 4000.,
    zcnt = 3000.,
 /
 &surface_flux
    ishflux = 0,
    tdelt = 12.,
    icoast = 30,
    cdh = 7.2e-3,
    irand = 0,
 /
 &streamfunction
    istrfcn = 0,
    s_repeat = 0,
    s_ampl = 40.,
    s_znaught = 6000.,
    s_hwavel = 10000.,
    s_vwavel = 18000.,
    s_period = 1200.,
 /
 &atmos_heat_source
    ihsrc = 0,
    h_ampl = 0.075,
    h_radius_x = 3000.,
    h_radius_z = 3000.,
    h_center_z = 3000.,
    h_freq = 0.005,
    h_modes = 2,
 /
 &rotunno_seabreeze
    iseabreeze = 0,
    sb_ampl = 0.000175,
    sb_x0 = 1000.,
    sb_z0 = 1000.,
    sb_period = 1.0,
    sb_latitude = 60.,
    sb_linear = 1,
 /
 &cooling_zone
    icoolzone = {icz},
    cz_ampl = 4.5,
    cz_rightedge = 0.20,
    cz_depth = 2000.,
    cz_width = 4000.,
    cz_coolrate = 600.,
 /
 """

    script_directory = os.path.dirname(__file__)

    # Specify the file path where you want to save the text file
    file_path = os.path.join(script_directory, "data.txt")

    # Write the user data to a text file
    with open(file_path, 'w') as file:
        file.write(data)
    result_label.config(text=f"User data has been written to {file_path}.")

def run_model(plot_option):
    global result_label  # Access the global variable
    # Get the directory of the current script
    script_directory = os.path.dirname(__file__)

    # Assuming your model executable is named 'dtdm'
    model_executable = os.path.join(script_directory, './dtdm')

    # Assuming the data file is named 'data.txt'
    data_file = os.path.join(script_directory, 'data.txt')

    # Run the model with the created configuration file
    model_command = f"{model_executable} < {data_file}"

    # Menjalankan model menggunakan subprocess
    try:
        subprocess.run(model_command, shell=True, check=True)
        result_label.config(text=f"Model has been successfully run using {data_file}.")

        # Menambahkan konversi menggunakan xgrads
        ctl_file_path = 'coldpool.hires.ctl'
        nc_file_path = 'coldpool.hires.nc'

        grads_ctl = xgrads.open_mfdataset(ctl_file_path)
        grads_ctl.to_netcdf(nc_file_path)

        print(f"File {ctl_file_path} berhasil dikonversi menjadi {nc_file_path}.")

    except subprocess.CalledProcessError as e:
        result_label.config(text=f"Error running the model or cdo command: {e}")

    if plot_option:
        show_plot()

def show_plot():
    dFil = "coldpool.hires.nc"
    fFil = os.path.join(dFil)

    ds = xr.open_dataset(fFil)
    cds = ds["thp"][-1, :, 0, :].squeeze()
    X, Y = np.meshgrid(cds["lon"].values, cds["lev"].values)
    z = cds.values

    # Create a contour plot using Matplotlib
    plt.figure(figsize=(8, 6))
    contour_plot = plt.contourf(X, Y, z, levels=20, cmap='plasma')

    # Add colorbar
    cbar = plt.colorbar(contour_plot, label="Temperature")

    # Set labels and title
    plt.xlabel("Longitude")
    plt.ylabel("Level")
    plt.title("Contour Plot of Temperature")

    # Show the plot
    plt.show()

def info():
    about_msg="""
------------------------------------------------------------------------------
                   INTRODUCTION
------------------------------------------------------------------------------
The grid_run namelist specifies model dimensions, integration length, and plotting output interval

 nx - number of horizontal grid points - max NXM in storage.txt default is 101
 nz - number of vertical grid points - max NZM in storage.txt default is 84
 dx - horizontal grid spacing (m; default = 1000.)
 dz - vertical grid spacing (m; default = 250.)
 dt - time step (s; default = 1.0)
 timend - integration stop time (s; default = 7200.)
 plot - plotting interval (s; default = 300.)
   * if interval < 60 sec, GrADS will report time incorrectly

===============================================================================
Framework namelist specifies whether the model is compressible or anelastic

 ipressure - output pressure decompositions if 1 (default is 0)
    * can make model integrations much longer in some cases
 imoist - turns on moisture if 1 (default is 0)

===============================================================================

The cooling_zone namelist implements a lower tropospheric cooling source following Fovell and Tan (2000) or an impulsive block of cold air of specified dimensions

 icoolzone (1 = Fovell-Tan storm-adaptive cooling zone,         2 = impulsive cold block; default is 0)
-------------------------------------------------------------------------------
    """
    top = tk.Toplevel(width=200, height=500)
    top.title("Mountain Lee Wave Model Help/Info")
    
    txt = tk.Text(top, borderwidth=3, relief=tk.RIDGE)
    txt.config(font=("consolas", 12), wrap='word')

    scrollbar = tk.Scrollbar(top, command=txt.yview)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    txt.insert(tk.END, about_msg)
    txt.config(state=tk.DISABLED)
    txt.pack()
    txt['yscrollcommand'] = scrollbar.set
    
    button = tk.Button(top, text="Close", command=top.destroy)
    button.pack()
    
    
    
# Create the main application window
root = tk.Tk()
app = FortranCompilerGUI(root)
root.mainloop()
