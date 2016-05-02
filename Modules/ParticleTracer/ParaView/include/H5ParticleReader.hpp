
#ifdef _cplusplus
extern "C" {
#endif
#include <hdf5.h>
#ifdef _cplusplus
}
#endif

#include <vector>

#define MAX_NAME 256

using namespace std;

class H5ParticleReader
{
    private:
        static void extract1DFloatData(hid_t group, char * datasetName, vector<float> * vData );

        /* Dipatches data initially packed in vX as x, y, z, x, y, z ...
         * By moving y coordinates to vY vector and z coordinates in vZ vector
         */
        static void dispatch(vector<float> * vX, vector<float> * vY, vector<float> * vZ);

        static void fillTimeSteps(vector<float> * vX, vector<float> * vTimeStep);

    public:
        static int readHDF5(char * filename, vector<vector<float> > * vvTimeStep, vector< vector<float> > * vvX, vector< vector<float> > * vvY, vector< vector<float> > * vvZ);
};

//int main(int argc, char ** argv)
//{
    ////readAll(argc, agrv);

    //vector<vector<float> > ts;
    //vector<vector<float> > x;
    //vector<vector<float> > y;
    //vector<vector<float> > z;

    //readHDF5("/tmp/PT/points.h5", &ts, &x, &y, &z);

    //printf("%d %d %d %d %d %d\n", ts.size(), ts.front().size(), x.size(), x.front().size(), y.size(), z.size());

    //return 0;
//}
