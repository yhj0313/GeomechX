#!/usr/bin/env python
# VTK high order: https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/

# modified from: $PETSC_DIR/lib/petsc/bin/petsc_gen_xdmf.py
import h5py
import numpy as np
import os, sys

class Xdmf:
  def __init__(self, filename):
    self.filename = filename
    self.cellMap  = {1 : {1 : 'Polyvertex', 2 : 'Polyline'}, 2 : {3 : 'Triangle', 4 : 'Quadrilateral'}, 3 : {4 : 'Tetrahedron', 6: 'Wedge', 8 : 'Hexahedron'}}

    # py2/py3 compatibility, see https://github.com/h5py/h5py/issues/379
    if sys.version_info[0] < 3:
      self.typeMap = {'scalar' : 'Scalar', 'vector' : 'Vector', 'tensor' : 'Tensor6', 'matrix' : 'Matrix'}
      self.typeExt = {2 : {'vector' : ['x', 'y'], 'tensor' : ['xx', 'yy', 'xy']}, 3 : {'vector' : ['x', 'y', 'z'], 'tensor' : ['xx', 'yy', 'zz', 'xy', 'yz', 'xz'], 'scalar' : ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']}}
    else:
      self.typeMap = {b'scalar' : 'Scalar', b'vector' : 'Vector', b'tensor' : 'Tensor6', b'matrix' : 'Matrix'}
      self.typeExt = {2 : {b'vector' : ['x', 'y'], b'tensor' : ['xx', 'yy', 'xy'], b'scalar' : ['xx', 'yy', 'xy', 'zz']}, 3 : {b'vector' : ['x', 'y', 'z'], b'tensor' : ['xx', 'yy', 'zz', 'xy', 'yz', 'xz'], b'scalar' : ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']}}
    return

  def writeHeader(self, fp, hdfFilename):
    fp.write('''\
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [
<!ENTITY HeavyData "%s">
]>
''' % os.path.basename(hdfFilename))
    fp.write('\n<Xdmf>\n  <Domain Name="domain">\n')
    return

  def writeCells(self, fp, topologyPath, numCells, numCorners, cellsName = "cells"):
    fp.write('''\
    <DataItem Name="%s"
              ItemType="Uniform"
              Format="HDF"
              NumberType="Float" Precision="8"
              Dimensions="%d %d">
      &HeavyData;:/%s/cells
    </DataItem>
''' % (cellsName, numCells, numCorners, topologyPath))
    return

  def writeVertices(self, fp, geometryPath, numVertices, spaceDim):
    fp.write('''\
    <DataItem Name="vertices"
              Format="HDF"
              Dimensions="%d %d">
      &HeavyData;:/%s/vertices
    </DataItem>
    <!-- ============================================================ -->
''' % (numVertices, spaceDim, geometryPath))
    return

  def writeLocations(self, fp, numParticles, spaceDim):
    fp.write('''\
    <DataItem Name="particle_coordinates"
              Format="HDF"
              Dimensions="%d %d">
      &HeavyData;:/particles/coordinates
    </DataItem>
''' % (numParticles, spaceDim))
    return

  def writeTimeGridHeader(self, fp, time):
    fp.write('''\
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
      <Time TimeType="List">
        <DataItem Format="XML" NumberType="Float" Dimensions="%d">
          ''' % (len(time)))
    fp.write(' '.join([str(float(t)) for t in time]))
    fp.write('''
        </DataItem>
      </Time>
''')
    return

  #http://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology
  def writeHybridSpaceGridHeader(self, fp):
    fp.write('      <Grid Name="domain" GridType="Collection">\n')
    return

  def writeSpaceGridHeader(self, fp, numCells, numCorners, cellDim, spaceDim, cellsName = "cells"):
    fp.write('''\
      <Grid Name="domain" GridType="Uniform">
        <Topology
           TopologyType="%s"
           NumberOfElements="%d">
          <DataItem Reference="XML">
            /Xdmf/Domain/DataItem[@Name="%s"]
          </DataItem>
        </Topology>
        <Geometry GeometryType="%s">
          <DataItem Reference="XML">
            /Xdmf/Domain/DataItem[@Name="vertices"]
          </DataItem>
        </Geometry>
''' % (self.cellMap[cellDim][numCorners], numCells, cellsName, "XYZ" if spaceDim > 2 else "XY"))
    return

  def writeFieldSingle(self, fp, numSteps, timestep, spaceDim, name, f, domain):
    if len(f[1].shape) > 2:
      dof = f[1].shape[1]
      bs  = f[1].shape[2]
    elif len(f[1].shape) > 1:
      if numSteps > 1:
        dof = f[1].shape[1]
        bs  = 1
      else:
        dof = f[1].shape[0]
        bs  = f[1].shape[1]
    else:
      dof = f[1].shape[0]
      bs  = 1
    
    if bs > 5:
      fp.write('''\
        <Attribute
           Name="%s"
           Type="%s"
           Center="%s">
          <DataItem ItemType="HyperSlab"
        	    Dimensions="1 %d %d"
        	    Type="HyperSlab">
            <DataItem
               Dimensions="3 3"
               Format="XML">
              %d 0 0
              1 1 1
              1 %d %d
            </DataItem>
            <DataItem
               DataType="Float" Precision="8"
               Dimensions="%d %d %d"
               Format="HDF">
              &HeavyData;:%s
            </DataItem>
          </DataItem>
        </Attribute>
  ''' % (f[0], 'Tensor6', domain, dof, bs, timestep, dof, bs, numSteps, dof, bs, name))
    else: 
      fp.write('''\
        <Attribute
           Name="%s"
           Type="%s"
           Center="%s">
          <DataItem ItemType="HyperSlab"
        	    Dimensions="1 %d %d"
        	    Type="HyperSlab">
            <DataItem
               Dimensions="3 3"
               Format="XML">
              %d 0 0
              1 1 1
              1 %d %d
            </DataItem>
            <DataItem
               DataType="Float" Precision="8"
               Dimensions="%d %d %d"
               Format="HDF">
              &HeavyData;:%s
            </DataItem>
          </DataItem>
        </Attribute>
  ''' % (f[0], self.typeMap[f[1].attrs['vector_field_type']], domain, dof, bs, timestep, dof, bs, numSteps, dof, bs, name))
    return

  def writeFieldComponents(self, fp, numSteps, timestep, spaceDim, name, f, domain):
    vtype = f[1].attrs['vector_field_type']
    if len(f[1].shape) > 2:
      dof    = f[1].shape[1]
      bs     = f[1].shape[2]
      cdims  = '1 %d 1' % dof
      dims   = '%d %d %d' % (numSteps, dof, bs)
      stride = '1 1 1'
      size   = '1 %d 1' % dof
    else:
      dof    = f[1].shape[0]
      bs     = f[1].shape[1]
      cdims  = '%d 1' % dof
      dims   = '%d %d' % (dof, bs)
      stride = '1 1'
      size   = '%d 1' % dof
    for c in range(bs):
      ext = self.typeExt[spaceDim][vtype][c]
      if len(f[1].shape) > 2: start  = '%d 0 %d' % (timestep, c)
      else:                   start  = '0 %d' % c
      fp.write('''\
        <Attribute
           Name="%s"
           Type="Scalar"
           Center="%s">
          <DataItem ItemType="HyperSlab"
        	    Dimensions="%s"
        	    Type="HyperSlab">
            <DataItem
               Dimensions="3 %d"
               Format="XML">
              %s
              %s
              %s
            </DataItem>
            <DataItem
               DataType="Float" Precision="8"
               Dimensions="%s"
               Format="HDF">
              &HeavyData;:%s
            </DataItem>
          </DataItem>
        </Attribute>
''' % (f[0]+'_'+ext, domain, cdims, len(f[1].shape), start, stride, size, dims, name))
    return

  def writeField(self, fp, numSteps, timestep, cellDim, spaceDim, name, f, domain):
    if sys.version_info[0] < 3:
      ctypes = ['tensor', 'matrix']
    else:
      ctypes = [b'tensor', b'matrix']
    if spaceDim == 2 or cellDim != spaceDim: 
      if sys.version_info[0] < 3: ctypes.append('vector')
      else: ctypes.append(b'vector')
    if f[1].attrs['vector_field_type'] in ctypes:
      self.writeFieldComponents(fp, numSteps, timestep, spaceDim, name, f, domain)
    elif f[1].attrs['vector_field_type'] == b'scalar':
      if not(numSteps < 2 and timestep == 0):
        if len(f[1].shape) > 2:
          self.writeFieldComponents(fp, numSteps, timestep, spaceDim, name, f, domain)
        else:
          self.writeFieldSingle(fp, numSteps, timestep, spaceDim, name, f, domain)
      else:
        if len(f[1].shape) > 1:
          self.writeFieldComponents(fp, numSteps, timestep, spaceDim, name, f, domain)
        else:
          self.writeFieldSingle(fp, numSteps, timestep, spaceDim, name, f, domain)
    return

  def writeSpaceGridFooter(self, fp):
    fp.write('      </Grid>\n')
    return

  def writeParticleGridHeader(self, fp, numParticles, spaceDim):
    fp.write('''\
      <Grid Name="particle_domain" GridType="Uniform">
        <Topology TopologyType="Polyvertex" NodesPerElement="%d" />
        <Geometry GeometryType="%s">
          <DataItem Reference="XML">/Xdmf/Domain/DataItem[@Name="particle_coordinates"]</DataItem>
        </Geometry>
''' % (numParticles, "XYZ" if spaceDim > 2 else "XY"))
    return

  def writeParticleField(self, fp, fieldname, numParticles, numComp):
    fp.write('''\
    <Attribute Name="particles/%s">
      <DataItem Name="%s"
                Format="HDF"
                Dimensions="%d %d">
                &HeavyData;:/particle_fields/%s
      </DataItem>
    </Attribute>
''' % (fieldname, fieldname, numParticles, numComp, fieldname))
    return

  def writeTimeGridFooter(self, fp):
    fp.write('    </Grid>\n')
    return

  def writeFooter(self, fp):
    fp.write('  </Domain>\n</Xdmf>\n')
    return

  def write(self, hdfFilename, topologyPath, numCells, numCorners, cellDim, htopologyPath, numHCells, numHCorners, geometryPath, numVertices, spaceDim, time, vfields, cfields, numParticles, pfields):
    useTime = not (len(time) < 2 and time[0] == -1)
    with open(self.filename, 'w') as fp:
      self.writeHeader(fp, hdfFilename)
      # Field information
      self.writeCells(fp, topologyPath, numCells, numCorners)
      if numHCells:
        self.writeCells(fp, htopologyPath, numHCells, numHCorners, "hcells")
      self.writeVertices(fp, geometryPath, numVertices, spaceDim)
      if useTime: self.writeTimeGridHeader(fp, time)
      for t in range(len(time)):
        if numHCells:
          self.writeHybridSpaceGridHeader(fp)
          self.writeSpaceGridHeader(fp, numHCells, numHCorners, cellDim, spaceDim, "hcells")
          self.writeSpaceGridFooter(fp)
        self.writeSpaceGridHeader(fp, numCells, numCorners, cellDim, spaceDim)
        for vf in vfields: self.writeField(fp, len(time), t, cellDim, spaceDim, '/vertex_fields/'+vf[0], vf, 'Node')
        for cf in cfields: self.writeField(fp, len(time), t, cellDim, spaceDim, '/cell_fields/'+cf[0], cf, 'Cell')
        self.writeSpaceGridFooter(fp)
        if numHCells:
          self.writeSpaceGridFooter(fp)
      if useTime: self.writeTimeGridFooter(fp)
      if numParticles:
        self.writeLocations(fp, numParticles, spaceDim)
        if useTime: self.writeTimeGridHeader(fp, time)
        for t in range(len(time)):
          self.writeParticleGridHeader(fp, numParticles, spaceDim)
          for pf in pfields:
            self.writeParticleField(fp, pf[0], numParticles, int(pf[1].attrs['Nc']))
          self.writeSpaceGridFooter(fp)
        if useTime: self.writeTimeGridFooter(fp)
      self.writeFooter(fp)
    return

def generateXdmf(hdfFilename, xdmfFilename = None):
  if xdmfFilename is None:
    xdmfFilename = os.path.splitext(hdfFilename)[0] + '.xmf'
  # Read mesh
  h5          = h5py.File(hdfFilename, 'r')
  if 'viz' in h5 and 'geometry' in h5['viz']:
    geomPath  = 'viz/geometry'
    geom      = h5['viz']['geometry']
  else:
    geomPath  = 'geometry'
    geom      = h5['geometry']
  if 'viz' in h5 and 'topology' in h5['viz']:
    topoPath  = 'viz/topology'
    topo      = h5['viz']['topology']
  else:
    topoPath  = 'topology'
    topo      = h5['topology']
  if 'viz' in h5 and 'hybrid_topology' in h5['viz']:
    htopoPath = 'viz/hybrid_topology'
    htopo     = h5['viz']['hybrid_topology']
  else:
    htopoPath = None
    htopo     = None
  vertices    = geom['vertices']
  numVertices = vertices.shape[0]
  spaceDim    = vertices.shape[1]
  cells       = topo['cells']
  numCells    = cells.shape[0]
  numCorners  = cells.shape[1]
  cellDim     = topo['cells'].attrs['cell_dim']
  if htopo:
    hcells      = htopo['cells']
    numHCells   = hcells.shape[0]
    numHCorners = hcells.shape[1]
  else:
    numHCells   = 0
    numHCorners = 0
  if 'time' in h5:
    time      = np.array(h5['time']).flatten()
  else:
    time      = [-1]
  vfields     = []
  cfields     = []
  pfields     = []
  pfields     = []
  if 'vertex_fields' in h5: vfields = h5['vertex_fields'].items()
  if 'cell_fields' in h5: cfields = h5['cell_fields'].items()
  numParticles = 0
  if 'particles' in h5:
    numParticles = h5['particles']['coordinates'].shape[0]
  if 'particle_fields' in h5: pfields = h5['particle_fields'].items()

  # Write Xdmf
  Xdmf(xdmfFilename).write(hdfFilename, topoPath, numCells, numCorners, cellDim, htopoPath, numHCells, numHCorners, geomPath, numVertices, spaceDim, time, vfields, cfields, numParticles, pfields)
  h5.close()
  return

if __name__ == '__main__':
  for f in sys.argv[1:]:
    generateXdmf(f)
