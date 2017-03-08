import vtktools

vtk_writer = vtktools.VTK_XML_Serial_Unstructured()

''' Only plot cells '''

x = [-2.0,0.0,1.5]
y = [2.0,-2.5,3.0]
z = [0.0, 0.0,0.0]
R = [2.0,2.5,3.0]



''' Plot cells and forces '''

F_x = [10.0,10.0,20.0]
F_y = [15.0,-5.0,0.0]
F_z = [0.0,0.0,0.0]

vtk_writer.snapshot("cell_arrangements.vtu", x, y, z, radii = R, x_force = F_x, y_force = F_y, z_force = F_z )
vtk_writer.writePVD("cell_arrangements.pvd")

#vtk_writer.snapshot("cell_arrangements.vtu", x, y, z, radii = R)
#vtk_writer.writePVD("cell_arrangements.pvd")

