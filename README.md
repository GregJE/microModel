# microModel
Greg Evans - 2018

microModel:
Calculates and Represents the diaphragm deformation in a microphone, given a acoustic pressure and tension force. This excludes the electrostatic interaction. Deformation is imposed as zero on boundaries.

Created as part of a Bachelors Degree (Hons) in Mechanical Engineering.

The following script solves the governing second-order differential equation through both analytical and numerical solution methodologies.
Equation: -T(d^2f/dx^2+d^f/dy^2) = Pa; (x,y) ? Lx x Ly

INPUT:
dx: increment size in x
dy: increment size in y
Lx: length of diaphragm in x
Ly: length of diaphragm in y
T: tension force
Pa: acoustic pressure

OUTPUT:
1. 3D Surface Model of the deformed diaphragm
2. x-y Comparison of numerical and analytical solution methods*
3. x-y Averaged cumulative error between solution methods
4. 2D Contour Plot
    *Only accurate for a square mesh
