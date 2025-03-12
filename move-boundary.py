import numpy as np

fp1 = "mach-adapted-mesh.su2"

nodes = []
cells = []
tags  = []

fpcell = "cellsize_old.log"
with open(fpcell, 'r') as cellsizefile:
    lines = cellsizefile.readlines()

    for line in lines:
        line = line.split('=')
        y0_tmp = line[1]
        break

y0 = float(y0_tmp)		# Input mesh first layer height

with open(fp1, 'r') as su2file:
    lines = su2file.readlines()
    line1 = lines[0].strip()
    if line1.startswith("NDIME="):
        dim = int(line1.split('=')[1])

    nodeSec = False
    cellSec = False
    tagSec  = False
    
    for line in lines[1:]:

        line = line.strip()       

        if line.startswith("NPOIN="):
            num_nodes = int(line.split('=')[1])
            nodeSec = True
            continue
        elif nodeSec and line.startswith("NELEM="):
            nodeSec = False

        
        if nodeSec:
            node_data = line.split()
            if len(node_data) == dim:
                nodes.append(list(map(float, node_data)))

        if line.startswith("NELEM="):
            k = 1
            cellSec = True
            num_cells = int(line.split('=')[1])
            continue
        elif cellSec and line.startswith("NMARK="):
            cellSec = False

        if cellSec:
            l = 1
            cell_data = line.split()
            if cell_data:
                cells.append(list(map(int, cell_data)))

        if line.startswith("NMARK="):
            tagSec = True
            num_tags = int(line.split('=')[1])
            continue

        if tagSec:
            tag_data = line.split()
            if tag_data:
                tags.append(list(map(str, tag_data)))
                
triscells_prim = []          # Contains all the tri cells
quadcells_prim = []          # Contains all the quad cells

for cell in cells:
    if cell[0] == 5:
        triscells_prim.append(cell)
    else:
        quadcells_prim.append(cell)

triscells = []              # Contains all the tri cells without the cell type marker
quadcells = []              # Contains all the quad cells without the cell type marker

for triscell in triscells_prim:
    triscells.append(triscell[1:])

for quadcell in quadcells_prim:
    quadcells.append(quadcell[1:])

qnodeInd = []           # Nodes that are associated with quad but not on the interface

for quadcell in quadcells:
    for quadnode in quadcell:
        if any(quadnode in row for row in triscells) == False and quadnode not in qnodeInd:
            qnodeInd.append(int(quadnode))

intnodesInd = []        # Nodes on interface only
intnodes = []           # Same as the above, repititions removed

for quadcell in quadcells:
    for quadnode in quadcell:
        if any(quadnode in row for row in triscells) == True and quadnode not in qnodeInd:
            intnodesInd.append(quadnode)

for intnodesIndi in intnodesInd:
    if intnodesIndi not in intnodes:
        intnodes.append(int(intnodesIndi))

inttriscells = []       # Triangular cells on the interface

for triscell in triscells:
    for trisnode in triscell:
        if trisnode in intnodes and triscell not in inttriscells:
            inttriscells.append(triscell)
            
walltag = '1'       # Name the tag assigned to the wall
wallcells_prim = []
wallcells = []
wallstart = False

for tag in tags:
    if wallstart == False and tag == ["MARKER_TAG=", walltag]:
        wallstart = True
        continue
    elif wallstart == True and tag[0] == 'MARKER_TAG=':
        wallstart = False
        continue
    elif wallstart:
        wallcells_prim.append(tag)
        continue
    elif wallstart == False:
        continue

for wallcell in wallcells_prim:
    wallcells.append(wallcell[1:])
    
wallnodesInd = []          # Nodes on the wall only
wallnodes    = []          # Nodes on the wall only, repititions removed

for wallcell in wallcells[1:]:
    for wallnode in wallcell:
        if wallnode not in wallnodesInd:
            wallnodesInd.append(int(wallnode))

for wallnodesIndi in wallnodesInd:
    if wallnodesIndi not in wallnodes:
        wallnodes.append(wallnodesIndi)
        
nodepairsInd = []
nodepairs = []          # Node pairs corresponding to wall and the interface, [wall, interface]

for wallnode in wallnodes:
    
    compnode = 0
    dist     = 0
    distnode = 0
    
    x_wall = nodes[wallnode][0]; x_wall = float(x_wall)
    y_wall = nodes[wallnode][1]; y_wall = float(y_wall)
    wall_coord = np.array((x_wall, y_wall))

    for intnode in intnodes:
 
        x_int = nodes[intnode][0]; x_int = float(x_int)
        y_int = nodes[intnode][1]; y_int = float(y_int)
        int_coord = np.array((x_int, y_int))
        
        dist = np.linalg.norm(int_coord - wall_coord)

        if distnode == 0:
            distnode = dist
            compnode = intnode
            continue
        elif distnode != 0 and dist < distnode:
            distnode = dist
            compnode = intnode
        else:
            continue

    nodepairsInd.append([wallnode, compnode])  

for nodepairsIndi in nodepairsInd:
    if nodepairsIndi not in nodepairs:
        nodepairs.append(nodepairsIndi)
        
innernodes = []         # Nodes in the interior of the quad region, i.e., not on the wall nor at the interface

for qnodeIndi in qnodeInd:
    if qnodeIndi not in wallnodes:
        innernodes.append(qnodeIndi)

## Hyperbolic tangent ##

def funct(delta, b):
   return np.sinh(delta) - b*delta

def functdash(delta, b):
   return np.cosh(delta) - b

def newtonraph(delta, b):
   h = 1.0
   while abs(h) >= 1e-15:
      h = funct(delta, b)/ functdash(delta, b)
      delta = delta - h
   return delta

def hyptangrid(delmin1, delmin2, p1, p2, nninit):

   delmin1 = delmin1/np.linalg.norm(p2 - p1)
   delmin2 = delmin2/np.linalg.norm(p2 - p1)
   nn = nninit
   a = np.sqrt(delmin2/delmin1)
   b = 1.0/(nn*np.sqrt(delmin1*delmin2))
   delta0 = 20.0
   delta = newtonraph(delta0, b)
   #delta = analyticdelta(b)
   s = np.zeros(nn)
   s[nn-2] = 0.99
   s[nn-3] = 0.99 - 2.0*delmin2
   # while (s[-2] - s[-3]) > delmin2:
   #    s = np.zeros(nn)
   #    u = np.zeros(nn)
   #    x = np.zeros(nn)
   #    eps = np.arange(nn)
   #    u[:] = 0.5*(1.0 + np.tanh(delta*(eps[:]/nn - 0.5))/np.tanh(0.5*delta))
   #    s[:] = u[:]/(a + (1.0 - a)*u[:])
   #    nn = nn + 10  

   s = np.zeros(nn)
   u = np.zeros(nn)
   x = np.zeros(nn)
   eps = np.arange(nn)
   u[:] = 0.5*(1.0 + np.tanh(delta*(eps[:]/nn - 0.5))/np.tanh(0.5*delta))
   s[:] = u[:]/(a + (1.0 - a)*u[:])
   return s
   
triareas = []       # Area of the triangles touching the interface
trisizes = []       # Quantitative approximation of the triangle size = (Area)^0.5
trisizesq = []      # Squares of the sizes (to find rms)

for intricell in inttriscells:
    p1 = np.array(nodes[intricell[0]])
    p2 = np.array(nodes[intricell[1]])
    p3 = np.array(nodes[intricell[2]])
    area = 0.5*np.abs((p1[0]*(p2[1]-p3[1]))+(p2[0]*(p3[1]-p1[1]))+(p3[0]*(p1[1]-p2[1])))
    triareas.append(area)

for triarea in triareas:
    size = (triarea)**0.5
    trisizes.append(size)

for trisize in trisizes:
    sizesq = trisize**2
    trisizesq.append(sizesq)

lastlayer = (np.average(np.array(trisizesq)))**0.5

fp2 = "paraview_files/yplus.sol"
yplus = []     # Y_plus at every wall node

with open(fp2, 'r') as surffile:
    lines = surffile.readlines()
    
    for line in lines[6:-1]:
        line = line.strip()
        yplus.append(float(line))
yplus = np.array(yplus)

flowfactor = yplus/y0
yd = []
for ff in flowfactor:
    yd.append(1/flowfactor)
firstlayer = np.min(yd)

i = 0
xy_quad = []        # New points to be placed in between each pair of wall and interface nodes
num_layers = 20
delmin1 = 1e-5
delmin2 = lastlayer

for nodepair in nodepairs:
    p1 = np.array(nodes[nodepair[0]])
    p2 = np.array(nodes[nodepair[1]])
    s  = hyptangrid(delmin1, delmin2, p1, p2, num_layers)
    new_points = p1 + (p2 - p1) * s[:, None]
    
    xy_quad.append(new_points[1:].tolist())
    
innerInds = np.zeros((len(xy_quad), len(xy_quad[0])))       # Node indices between every pair of wall-interface nodes corresponding to the sequence in "innernodes"
i = 0

for nodepair in nodepairs:
    a = np.array(nodes[nodepair[0]])
    b = np.array(nodes[nodepair[1]])

    j = 0

    for innernode in innernodes:
        p = np.array(nodes[innernode])

        cp = np.cross(p - a, b - a)
        cp = cp/(np.linalg.norm(p-a)*np.linalg.norm(b-a))
        #print(cp)

        if np.abs(cp) <= 0.01 and np.linalg.norm(p-a) < np.linalg.norm(b-a):
            innerInds[i][j] = innernode
            j += 1
        else:
            continue

    i += 1
    
k = 0
insortnodes = np.zeros(np.shape(innerInds))         # Same as "innerInds", but every row is re-arranged such that nodes start from wall and end at interface

for innerInd in innerInds:
    a = np.array(nodes[nodepairs[k][0]])
    distances = {i: np.linalg.norm(np.array(nodes[int(i)]) - np.array(a)) for i in innerInd}
    sortedInd = sorted(innerInd, key=lambda i: distances[i])
    insortnodes[k] = sortedInd
    k += 1

i = 0
for insortnode in insortnodes:
    j = 0
    for sortnode in insortnode:
        nodes[int(sortnode)] = xy_quad[i][j]
        j += 1
    i += 1

linesMod = []
linesMod.append("NDIME= " + str(dim))
linesMod.append("NPOIN= " + str(len(nodes)))
for node in nodes:
    linesMod.append(node)
linesMod.append("NELEM= " + str(len(cells)))
for cell in cells:
    linesMod.append(cell)
linesMod.append("NMARK= " + str(num_tags))
for tag in tags:
    linesMod.append(tag)

with open("mach-adapted-mesh-mod.su2", 'w') as file:
    for line in linesMod:
        if isinstance(line, list):  
            file.write(" ".join(map(str, line)) + '\n')
        else:  
            file.write(str(line) + '\n')
            
with open("cellsize.log", 'w') as sizefile:
    sizefile.write("firstlayer height = " + str(firstlayer) + '\n')
    sizefile.write("lastlayer height = " + str(lastlayer) + '\n')
