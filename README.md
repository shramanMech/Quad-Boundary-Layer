# Quad-Boundary-Layer
This work revolves around adjusting the structured boundary layer of a hybrid mesh to obtain better wall solution predictions in CFD analyses.

# Introduction
A hybrid mesh (unstructured mesh with a structured boundary layer) can be very useful for CFD simulations, especially if the flow is at hypersonic speeds, as it has been shown to improve wall heat flux predictions significantly better than an all-unstructured mesh [(Nastac et al.)](https://doi.org/10.2514/6.2022-0111), although its reliability against a completely structured mesh is still under study. Here, we use [NASA's Refine](https://github.com/nasa/refine?tab=readme-ov-file) to perform the unstructured mesh adaptation. While the Refine code adapts the unstructured part of the mesh (tets, in case of 3D, and triangles, in case of 2D) very effectively [(Prof. Balan's Git repo)](https://github.com/aravind-balan/Mesh-Adaptation/tree/main), it keeps the structured boundary part of the mesh fixed. Work is going on for quad adaptation strategies, but for the time being, we have worked on a makeshift procedure to adjust the structured boundary layer based on the flow solution. 
For more information on how metric field-based mesh adaptation is carried out, please visit [Prof. Balan's Git repo](https://github.com/aravind-balan/Mesh-Adaptation/tree/main).

Here, the Python code *move-boundary.py* takes in a hybrid mesh in *.su2* format and adjusts the height of the layers parallel to the streamwise direction. The adjustment for the first layer height is made based on the y+ values obtained after every simulation (as the entire mesh adaptation workflow involves multiple simulations). The final layer height is obtained by finding out the average size of the triangles/tets at the interface of the quad boundary and the unstructured region to ensure a smooth transition in the sizes of the cells. The sizes of the layers between the first and the last layer transition smoothly using a [hyperbolic tangent function](https://www.cfd-online.com/Wiki/Structured_mesh_generation).

**Limitations:** The code only adjusts the quad layers by moving them. It does not add layers to ensure an optimal number of layers. This is a WIP and will be added soon. 

# Test Cases 

**2D HEG Cylinder**
For this case, a RANS simulation was performed on a cylinder using METACOMP's CFD++ solver. The original experimental study can be found here: [HEG Cylinder Experiment](https://arc.aiaa.org/doi/10.2514/6.2003-4252). 
The simulation was performed considering a viscous, laminar and chemical and thermal non-equilibrium flow. A similar simulation was performed using SU2-NEMO that can be found [here](https://www.mdpi.com/2226-4310/8/7/193). 

The initial mesh had *20577 cells (17737 tri and 2840 quad)*. In the final adaptation (9th run), the number of triangles was around 74000, and the number of quads was fixed at 2840. 

![00_noquad_mesh](https://github.com/user-attachments/assets/6132d02c-5097-4208-ab4f-a7d347c3fb98)
Initial mesh. 

![00_noquad_mach](https://github.com/user-attachments/assets/69336ace-e205-4d6d-bf71-f25847487527)
Mach contour on the initial mesh. 

![09_noquad_mesh](https://github.com/user-attachments/assets/5e858afe-52b3-45b9-9a38-de6120495c41)
Final mesh

![09_noquad_mach](https://github.com/user-attachments/assets/98b8f7c4-798e-49b1-83a5-ba0cc04b67e3)
Mach contour on the final mesh

![09_quad_closeup](https://github.com/user-attachments/assets/94e5e4f5-cc1b-45d6-be4f-b889b1a8f0fa)
Close-up at the boundary and shock

