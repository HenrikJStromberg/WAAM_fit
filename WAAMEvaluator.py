import gmsh
import numpy as np
from scipy.optimize import root
import os
import sys
import random
import tempfile
import platform
import subprocess
from stl import mesh
from typing import Tuple

gmsh.initialize()


def evaluateSpheres(input, output, triangulationSizing=0.0):
    """Evaluate all inner and outer spheres on given shape and save to output

    Args:
        input (string): Path to input shape
        output (string): Path to output shape
        triangulationSizing (float, optional): controls size of triangulation. Defaults to 0.0 for auto-sizing.
    """
    nc, inz, N = getTriangulation(input, triangulationSizing)
    cnts = nc[inz].mean(axis=1)  # center_points

    r = getSphereRadii(nc, inz, N)  # radii
    gradient = np.zeros_like(r)
    for i in range(gradient.shape[0]):
        neighbours = (np.isin(inz, inz[i]).sum(axis=1) == 2).nonzero()[0]
        neighbours = neighbours[r[neighbours] != -1]  # neglect invalid radii
        gradient[i] = np.linalg.norm(
            (r[i] - r[neighbours]) / np.linalg.norm(cnts[neighbours] - cnts[i], axis=1))
    r_scaled = r/r.max()
    grad_scaled = gradient - gradient.min()

    btm_95_percent = (grad_scaled < grad_scaled.max() * 0.95)
    grad_scaled[grad_scaled >= grad_scaled.max() * 0.95] = grad_scaled[btm_95_percent.nonzero()[0][grad_scaled[btm_95_percent].argmax()]]
    grad_scaled = grad_scaled / grad_scaled.max()

    # Save and show data in gmsh GUI:
    _, elementTags, __ = gmsh.model.mesh.getElements(2)
    views = []
    views.append(__add_as_view_to_gmsh__(elementTags[0].tolist(), r_scaled.tolist(), "Sphere Radii")) # type: ignore
    views.append(__add_as_view_to_gmsh__(elementTags[0].tolist(), grad_scaled.tolist(), "Radii Gradients")) # type: ignore

    if not os.path.exists(os.path.dirname(output)):
        os.mkdir(os.path.dirname(output))

    for v in views:
        # Set a green to blue color map as the ColorTable = {Green, Red} option is not yet available in API
        gmsh.view.option.set_number(v, "ColormapNumber", 17)
        gmsh.view.option.set_number(v, "ColormapSwap", 1)
        
        # Save views:
        file_names = ["radii_scaled", "gradient_scaled"]
        gmsh.view.write(v, output + file_names.pop() + ".msh")

    # Show last view by default:
    gmsh.view.option.set_number(views[-1], "Visible", 1)

    gmsh.fltk.initialize()
    while gmsh.fltk.is_available():
        gmsh.fltk.wait()
    gmsh.finalize()


def getTriangulation(input: str, triangulationSizing=0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create triangulation mesh on input file and return mesh

    The mesh will be returned as BREP (boundary representation) with the node coordinates, the inzidenz_matrix (which gives the relation between nodes and edges) and ?
    Args:
        input (string): file name with path (either .step/.stp or .stl are supported)
        triangulationSizing (float, optional): Controls mesh-sizing. Defaults to 0.0 for auto-sizing.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: node_coordinates, inzidenz_matrix, vector?
    """
    file_extension = input.split('.')[-1]
    if file_extension in ['stp', 'step']:
        gmsh.model.occ.importShapes(input)
        gmsh.model.occ.synchronize()
        if triangulationSizing != 0:
            gmsh.model.mesh.setSize(
                # Pass mesh sizing trough to points (entities with dimension 0)
                gmsh.model.occ.getEntities(0), triangulationSizing)
        gmsh.model.mesh.generate(2)
        nc, inz, N = __MshFromGmsh__()
    elif file_extension == 'stl':
        meshObj = mesh.Mesh.from_file(input)
        minEdLen = np.stack([np.linalg.norm(meshObj.v0 - meshObj.v1, axis=1),
                             np.linalg.norm(meshObj.v0 - meshObj.v2, axis=1),
                             np.linalg.norm(meshObj.v1 - meshObj.v2, axis=1)]).min()
        ncRaw = np.vstack([meshObj.v0, meshObj.v1, meshObj.v2])
        inz = np.arange(ncRaw.shape[0]).reshape((-1, 3))
        pntIsRenamed = np.zeros(ncRaw.shape[0], dtype=bool)
        nc = np.zeros([0, 3])
        N = -meshObj.normals
        N = N / np.tile(np.linalg.norm(N, axis=1), (3, 1)).T

        for i in range(ncRaw.shape[0]):
            if pntIsRenamed[i]:
                continue
            dists = np.linalg.norm(ncRaw - ncRaw[i], axis=1)
            identicalPnts = (dists < minEdLen / 2).nonzero()[0]
            pntIsRenamed[identicalPnts] = True
            pntIsRenamed[i] = False
            inz[np.isin(inz, identicalPnts)] = nc.shape[0]
            nc = np.vstack([nc, ncRaw[i]])
    else:
        raise ValueError('File format is not supported')
    return nc, inz, N


def getSphereRadii(nc: np.ndarray, inz: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Computes all radii of all spheres, indicated by the centers of the coordinate nodes. Invalid radii have the value -1.

    Args:
        nc (np.ndarray): Array with node coordinates
        inz (np.ndarray): Inzidenz matrix or vector?
        N (np.ndarray): Vector with unknown meaning

    Returns:
        np.ndarray: Vector with computed radii. Default values are -1. Values will only be changed, if a radius was sucessfully computed.
    """
    cnts = nc[inz].mean(axis=1)  # center points

    r = -np.ones(inz.shape[0])
    for i in range(inz.shape[0]):
        def f(r): return __evalRadius__(i, cnts, N, r)
        start = 2
        for j in range(10):
            start = start / 2
            solObj = root(f, start)
            if solObj['success']:
                r[i] = solObj['x']
                break
    return r


def plotCurrentGeo():
    fp = os.path.join(tempfile.gettempdir(), '0') + '.step'
    while os.path.exists(fp):
        fp = os.path.join(tempfile.gettempdir(), str(
            random.randint(0, 10 ** 6))) + '.step'
    gmsh.write(fp)
    subprocess.Popen(
        [sys.executable, os.path.join(os.path.abspath(
            os.getcwd()), 'GmshPlotter.py'), fp],
        shell=False)


    """Computes the distance from the point at index (current base point) to the center point other than index closest to a Heuver's sphere of radius r to the base point.

    Args:
        index (integer): index of current base point
        cnts (np.ndarray): nx3 matrix of center points 
        N (np.ndarray): nx3 matrix of normal vectors
        r (float): radius

    Returns:
        float: distance
    """
    basePnt = cnts[index]
    distances = np.linalg.norm((basePnt - cnts) * N, axis=1)
    distances[index] = np.inf
    return np.abs(distances-r).min()


def evaluateIslands(input, output, N=np.array([0, 0, 1]), stepSize=0):
    nc, inz, _ = __getMsh__(input)
    bb = np.reshape(gmsh.model.getBoundingBox(-1, -1), (2, 3))
    cnt = bb.mean(axis=0)
    N = N / np.linalg.norm(N)
    X0 = N * np.inner(cnt - bb[0], N)
    pathLen = np.linalg.norm(bb[1] - bb[0])
    if stepSize <= 0:
        ed = np.vstack([inz[:, [0, 1]], inz[:, [0, 2]], inz[:, [1, 2]]]).T
        stepSize = np.linalg.norm(nc[ed[0]] - nc[ed[1]], axis=1).min()
    steps = int(np.round(pathLen / stepSize))
    stepSize = pathLen/steps
    for i in range(steps):
        X = X0 + N * i * stepSize
        slice = __getSlice(nc, inz, X, N)

    gmsh.model.remove()


def __getSlice(nc, inz, X, N):
    ncOut = np.array([0, 3])
    edOut = np.array([0, 2], dtype=int)

    for i in range(inz.shape[0]):
        trgNodeDir = np.sign(np.inner(X - nc[inz[i]], N))
        # Is triangle is cut by plane X, N?
        if not (trgNodeDir == trgNodeDir[0]).all():
            A = nc[inz[i, 0]]
            B = nc[inz[i, 1]]
            C = nc[inz[i, 2]]

            if trgNodeDir[0] == trgNodeDir[1]:
                print('')  # ToDo: complete
            elif trgNodeDir[0] == trgNodeDir[2]:
                print('')  # ToDo: complete

            elif trgNodeDir[1] == trgNodeDir[2]:
                print('')  # ToDo: complete
    print('')


def __getMsh__(input):
    if input[-3:] == 'stp' or input[-4:] == 'step':
        gmsh.model.occ.importShapes(input)
        gmsh.model.occ.synchronize()
    elif input[-3:] == 'stl':
        gmsh.open(input)
    else:
        raise ValueError('unknown format')
    return __MshFromGmsh__()


def __MshFromGmsh__():
    """Do unknown operation

    Returns:
        tuple: nodeCoordinates(np.ndarray), elementNodeTags(np.ndarray), ?(np.ndarray)
    """
    nodeCoords = np.array([], dtype=float)
    gmsh.model.mesh.renumber_nodes()
    for i in range(0, 3):
        rnodeCoords = gmsh.model.mesh.getNodes(i)[1]
        nodeCoords = np.append(nodeCoords, rnodeCoords)
    nc = nodeCoords.reshape(-1, 3)
    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)
    elemTags = elemTags[0].astype(int)  # type: ignore
    inz = elemNodeTags[0].astype(int).reshape(  # type: ignore
        elemTags.shape[0], -1) - 1

    C = np.mean(nc[inz], axis=1)
    faceIDs = np.zeros(inz.shape[0], dtype=int)
    for entity in gmsh.model.getEntities(2):
        ID = entity[1]
        elemTagsOnFace = gmsh.model.mesh.getElements(2, ID)[1][0]
        for tag in elemTagsOnFace:
            faceIDs[elemTags == tag] = ID

    N = np.cross(nc[inz[:, 1]] - nc[inz[:, 0]], nc[inz[:, 2]] - nc[inz[:, 0]])
    for i in range(inz.shape[0]):
        gmsh.model.mesh.getNodes(-1, -1, returnParametricCoord=True)
        para = gmsh.model.getParametrization(2, faceIDs[i], C[i])
        N[i] = gmsh.model.getNormal(faceIDs[i], para)
        N[i] = N[i] / np.linalg.norm(N[i])

    return nc, inz, N


def __exportToOpenSCAD__(msh, outPath, elemNames=None, colors=None):
    if colors is list:
        colors = np.array(colors)
    if colors is not None and colors.ndim == 1:
        colors = np.expand_dims(colors, 1)

    elemTypeShortNames = {2: 'Triangle', 3: 'Quad',
                          4: 'Tetrahedron', 5: 'Hexahedron', 6: 'Prism', 7: 'Pyramid'}
    elemColors = {2: 'red', 3: 'blue', 4: 'blue',
                  5: 'red', 6: 'green', 7: 'yellow'}

    fileContent = ['//Mesh exported using MeshTools.exportToOpenSCAD\n']

    nc = msh['nc']

    # write point coordinates
    fileContent += ['Points = [\n']
    for i in range(nc.shape[0]):
        fileContent += [str(nc[i].tolist()) + ', // ' + str(i) + '\n']
    fileContent += ['  ];\n']

    # write elements
    e = 0
    for etIndex in range(msh['elemTypes'].shape[0]):
        et = msh['elemTypes'][etIndex]
        etFaces = [np.array([0, 1, 2], dtype=int)]
        inz = msh['inz'][etIndex]
        for i in range(inz.shape[0]):
            if elemNames is None:
                elemName = elemTypeShortNames[et] + str(i)
            else:
                elemName = str(elemNames[et][i])
            fileContent += ['\n' + elemName + ' = [\n']
            cnt = msh['nc'][inz[0]].mean(axis=0)
            for face in etFaces:
                N = np.cross(nc[inz[i, face[1]]] - nc[inz[i, face[0]]],
                             nc[inz[i, face[2]]] - nc[inz[i, face[1]]])
                if np.inner(np.mean(nc[inz[i, face]], axis=0) - cnt, N) > 0:
                    fileContent += [str(inz[i, face].tolist()) + ',']
                else:
                    fileContent += [str(inz[i, np.flip(face)].tolist()) + ',']
            fileContent[-1] = fileContent[-1][:-1]
            fileContent += ['  ];\n']

            fileContent += ['//[' + str(inz[i].tolist()) + ']\n']
            if colors is None:
                fileContent += ['color("' + elemColors[et] +
                                '") polyhedron( Points, ' + elemName + ' );\n']
            else:
                if colors.shape[1] == 3:
                    colorEntry = '[' + str(colors[e, 0]) + ', ' + \
                        str(colors[e, 0]) + ', ' + str(1 - colors[e]) + ']'
                else:
                    if colors[e] >= 0:
                        colorEntry = '[' + str(colors[e, 0]) + \
                            ', 0., ' + str(1 - colors[e, 0]) + ']'
                    else:
                        colorEntry = '[1., 1., 1.]'
                fileContent += ['color(' + colorEntry +
                                ') polyhedron( Points, ' + elemName + ' );\n']
            e = e + 1

    fileContent += ['LineWidth = 0.03;\n']

    fileContent += ['module line(start, end, thickness) {\n']
    fileContent += ['    color("black") hull() {\n']
    fileContent += ['        translate(start) sphere(thickness);\n']
    fileContent += ['        translate(end) sphere(thickness);\n']
    fileContent += ['    }\n']
    fileContent += ['}\n']

    if not outPath[-5:] == '.scad':
        outPath += '.scad'
    fp = open(outPath, 'w')
    fp.writelines(fileContent)
    fp.close()


def plotSolid(nc, inz, value, autoLaunch=True):
    scadPath = os.path.join(tempfile.gettempdir(), 'out.scad')
    value = value - value.min()
    value = value/value.max()
    clormap = [''] * inz.shape[0]
    for i in range(inz.shape[0]):
        clormap[i] = '[' + str(value[i]) + ', 0., ' + str(1 - value[i]) + ']'

    __exportToOpenSCAD__(
        {'nc': nc, 'inz': [inz], 'elemTypes': np.array([2])}, scadPath, value)
    if autoLaunch:
        if platform.system() == 'Darwin':
            subprocess.call(('open', scadPath))
        elif platform.system() == 'Windows':
            os.startfile(scadPath)  # type: ignore
        else:
            subprocess.call(('xdg-open', scadPath))


def __add_as_view_to_gmsh__(tags, data, view_name):
    """
    """
    view = gmsh.view.add(view_name)
    gmsh.view.add_homogeneous_model_data(view, 0, "", "ElementData", tags=tags, data=data, numComponents=1)
    gmsh.view.option.set_number(view, "Visible", 0)
    return view