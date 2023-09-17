# Fully discrete Lyapunov consistent discretizations of any order for parabolic reaction-diffusion equations with r species

> This README file

The fully-discrete Lyapunov consistent discretizations are designed for systems of parabolic reaction-diffusion equations. These discretizations preserve the stability of the equilibrium points of the continuous partial differential equations. The spatial discretization is based on summation-by-parts operators and simultaneous approximation terms for unstructured grids, while the temporal integration relies on relaxation Runge-Kutta methods. The discretizations are arbitrarily high-order accurate in space and time and can be used with h- and p-refinements. The algorithms have been successfully applied to epidemiology compartmental models for disease transmission and an oncolytic M1 virotherapy model. The source codes to simulate the susceptible–infected (SI) PDE model, susceptible-exposed-infectious (SEI) PDE model, and the Oncolytic M1 virotherapy PDE model reported in the manuscript can be found in this repository.

In the following, we guide the command line options for utilizing the SSDC library integrated with PETSc to execute the codes. Users can leverage these options to customize various aspects such as mesh handling, visualization, time stepping, adaptivity, and monitoring.

# Code Execution

To run the codes, follow the steps below:

##  Run the Code
* Execute the command ```make run ``` to run all the codes (Susceptible–infected (SI) PDE model, Susceptible-exposed-infectious (SEI) PDE model, and Oncolytic M1 virotherapy PDE model).
* You can pass additional options to customize the behavior of the code. Refer to the options listed in the available customization options section.

## Customization Options

The code provides various customization options that can be passed as command-line arguments. Below are the available options categorized according to their functionality:

1. Generic Options

`-ssdc_name`: assign a name to the SSDC object.

`-ssdc_view`: view SSDC context.

`-ssdc_deg`: specify the degree of the solution, OptionsAlias ```-deg```.

2.  Builtin box mesh generator

`-ssdc_mesh_box_dim`: number of space dimensions, OptionsAlias ```-dim```.

`-ssdc_mesh_box_size`: number of cells per direction, OptionsAlias ```-nel```.

`-ssdc_mesh_box_wrap`: periodicity per direction, OptionsAlias ```-W```.

`-ssdc_mesh_box_bbox`: bounding box `xmin,xmax,ymin,ymax,zmin,zmax`, OptionsAlias ```-B```.


3. Reading external meshes

`-ssdc_mesh_read`: read mesh from file, OptionsAlias ```-read```.


4. Uniform refinement

`-ssdc_refine_uniform`: levels of uniform refinement, OptionsAlias ```-refine```.


5.  Time Stepping

  * TS Solver

     `-ts_dt`: specify the time step size.

     `-ts_init_time`: specify the initial time.

     `-ts_max_time`: specify the maximum simulation time.

     `-ts_max_steps`: specify the maximum number of time steps.

     `-ssdc_load_step`: specify the load step.

     `-ssdc_save_step`: specify the save step.

     `-ssdc_load_initial`: load the initial conditions.

     `-ssdc_save_solution`: save the solution.

     `-ssdc_view_initial`: view the initial conditions.

     `-ssdc_view_solution`: view the solution.


  * Runge-Kutta Methods

     `-ts_rk_type`: choose the Runge-Kutta scheme.

   * Adaptivity

     `-ts_atol`: set the absolute tolerance for adaptivity.

     `-ts_rtol`: set the relative tolerance for adaptivity.

     `-ts_adapt_type`: set the type of adaptivity.

     `-ts_adapt_dt_min`: set the minimum time step size for adaptivity.

     `-ts_adapt_dt_max`: set the maximum time step size for adaptivity.

6. Monitoring

`-ts_monitor`: enable time stepping monitoring.

`-ts_adapt_monitor`: enable adaptivity monitoring.

`-ssdc_monitor_info`: monitor general information.

`-ssdc_monitor_step`: monitor step information.

`-ssdc_monitor_vtk`: monitor VTK output.

`-ssdc_monitor_dt`: monitor time step size.


7. Specific Options

`-example`: select the value of the epidemiological parameter R0 from five different values for the SI PDE model.

`-case`: select the equilibrium points case for the oncolytic M1 virotherapy PDE model. 

`-eq`: print the equilibrium value.

`-max_cases_monitor`: print temporal evolution of the maximum norm of the solution.

`-rrk_monitor`: print evolution of the Lyapunov functional and its time derivative.

