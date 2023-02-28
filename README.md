# DSMC (Direct Simulation Monte Carlo)

# Example

Heat diffusion through Argon gas at STP when the bottom boundary is hot:
* There are 10^5 simulated particles, each representing ~270 million Argon atoms
* The side length of the cubic domain is 100 microns
* The computational grid is 2x2
* Reflection off of the domain boundaries is diffuse

Initial snapshot  |  . |  .| .
:-------------------------:|:-------------------------:|:-----------:|---------------:
![image](https://user-images.githubusercontent.com/49919175/221903145-bac8047f-21c9-45de-88ff-ddcd78b39015.png)  |  ![image](https://user-images.githubusercontent.com/49919175/221903182-cca674c3-abd5-4006-af05-1fca7cdef571.png) | ![image](https://user-images.githubusercontent.com/49919175/221903233-a3799aa3-87a3-49f3-87b7-924b0ec8faf9.png) | ![image](https://user-images.githubusercontent.com/49919175/221903292-aac2b597-c7b0-4e08-b09b-366759c60e35.png)

# Overview

This project aims to use Direct Simulation Monte Carlo (DSMC) methods for modeling (dilute) gas flow.

The goal is to later extend this project to support simulation of dense gases / fluids.

See https://enchanting-turner-020.notion.site/DSMC-Readme-fe1d611d1d7d49e09fda3d6477bea094 for mathematical background and more information.
