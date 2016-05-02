# Particle Tracer

## Online Particle Tracer

(Example using the phantom dataset)

* Launch ParaView
* Load the dataset to use with the ParticleTracer
* Create a point source with:
  * Point coordinates: (36.621, 4.13897, 11.3298)
  * Number of points: >200
  * Radius: 2.0
* This will create a point source at the entry of the Phantom
* Select the initial data in the pipeline explorer
* Create a Particle Tracer filter
  * Select the initial data as input
  * Select the point source as seed source
* Change the Force reinjection every time steps (ex: every 10 timesteps)

* To save the particle data, add a Programmable filter with the following code as output for the Particle Tracer:

```
import os
import h5py

# Get input data
pdi = self.GetPolyDataInput()
pts = pdi.GetPointData()
nbPoints = pdi.GetNumberOfPoints()

# Get the ids of the particles
# They are deleted as they disappear
# from the geometry
ptids = pts.GetArray("ParticleId")

# Create a directory to store the
# particles for JEMRIS
if(not os.path.exists("/tmp/PT")):
    os.mkdir("/tmp/PT")
    
# Get the current timestep
ts = pdi.GetInformation().Get(output.DATA_TIME_STEP())

# Open the hdf5 file for storage
h5fl = h5py.File("/tmp/PT/points.h5", "a")

# Iterate over the particles and
# store them in .dat files
for i in range(nbPoints):
    coord = pdi.GetPoint(i)
    x, y, z = coord[:3]

    # The file will be structured as follows:
    # /<particle_id>/starttime: starting timestep
    # /<particle_id>/coords: successive coordinates of the particle
   
    grp = None
    
    # Try to open the dataset for the current particle
    # If the particle entry does not exist, we create it
    # and add the starttime initial data
    try:
        grp = h5fl[str(i)]
    except KeyError, e:
        grp = h5fl.create_group(str(i))

        dset = grp.create_dataset("starttime", (1,), 'f')
        dset[0] = ts

	 # Append to or create the coords field
    try:
        dset = grp["coords"]
        dset.resize((dset.len() + 3,))
    except KeyError, e:
        dset = grp.create_dataset("coords", (3,), maxshape=(None,), chunks=(500,))

    dset[dset.len() - 1] = x
    dset[dset.len() - 2] = y
    dset[dset.len() - 3] = z

h5fl.close()
```

* Click the play button to launch the particle tracer
* Once all the timesteps have been played, the data can be read in `/tmp/PT/points.h5`


## Offline particle tracer

* First, Build the `particletracer_pv` application, to do so:
  * Build AngioTK by adding `-DBUILD_MODULE_ParticleTracer=ON` and set an installation prefix with `-DCMAKE_INSTALL_PREFIX=...`;
  * Perform a `make install`;
  * The application should be in `<install_dir>/Modules/ParticleTracer/bin`.
* To use the application, you can specify the following options:
  * `--input.pointsources.filename <file.csv>`: a csv file containing the description of point sources to use. The format is: `position_x,position_y,position_z,radius,nbPoints`. An example is given for the Phantom dataset in the test directory;   
  * `--input.data.filename <file.case>`: The dataset on which you want to perform particle tracing;
  * `--output.vtk.enable 1`: By default, the particle tracer will only output an HDF5 file with particle data. With this option, you can additionnally add an export to vtk data to ensure the correctness of the results.

  
# Import HDF5 particles

* To import the particles generated in the HDF5 format, you can use the `H5ParticleReader` installed at the same time as the particle tracer. The usage is quite simple:
  * Link your application with the library
  * Use the provided header to access relevant functions
  * Import the data using the following static function (An example is given in the test directory: test_reader):
  
```
# time data
vector<vector<float> > ts;
# x positions
vector<vector<float> > x;
# y positions
vector<vector<float> > y;
# z positions
vector<vector<float> > z;

# Import data from hdf5 file into the vectors
H5ParticleReader::readHDF5("/tmp/PT/points.h5", &ts, &x, &y, &z);
```