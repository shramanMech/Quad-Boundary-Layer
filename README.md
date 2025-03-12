# Quad-Boundary-Layer
This work revolves around adjusting the structured boundary layer of a hybrid mesh (unstructured mesh with a structured boundary layer) to obtain better wall solution predictions while performing a CFD analysis.

# Introduction
A hybrid mesh can be very useful for CFD simulations, especially if the flow is at hypersonic speeds, as it has been shown to improve wall heat flux predictions significantly better than an all-unstructured mesh [(Nastac et al.)](https://doi.org/10.2514/6.2022-0111), although its reliability against a completely structured mesh is still under study. Here, we use [NASA's Refine](https://github.com/nasa/refine?tab=readme-ov-file) to perform the unstructured mesh adaptation. While the Refine code adapts the unstructured part of the mesh (tets, in case of 3D, and triangles, in case of 2D) very effectively [(Prof. Balan's Git repo)](https://github.com/aravind-balan/Mesh-Adaptation/tree/main), it keeps the structured boundary part of the mesh fixed. Work is going on for quad adaptation strategies, but for the time being, we have worked on a makeshift procedure to adjust the structured boundary layer based on the flow solution. 
For more information on how metric field-based mesh adaptation is carried out, please visit [Prof. Balan's Git repo](https://github.com/aravind-balan/Mesh-Adaptation/tree/main).

Here, the Python code *move-boundary.py* takes in a hybrid mesh in *.su2* format and adjusts the height of the layers parallel to the streamwise direction. The adjustment for the first layer height is made based on the y+ values obtained after every simulation (as the entire mesh adaptation workflow involves multiple simulations). The final layer height is obtained by finding out the average size of the triangles/tets at the interface of the quad boundary and the unstructured region to ensure a smooth transition in the sizes of the cells. The sizes of the layers between the first and the last layer transition smoothly using a [hyperbolic tangent function](https://www.cfd-online.com/Wiki/Structured_mesh_generation).

**Limitations:** The code only adjusts the quad layers by moving them. It does not add layers to ensure an optimal number of layers. This is a WIP and will be added soon. 

*test cases will be added soon*
