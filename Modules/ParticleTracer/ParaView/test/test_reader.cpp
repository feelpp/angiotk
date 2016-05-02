
#include "H5ParticleReader.hpp"

int main(int argc, char ** argv)
{
    vector<vector<float> > ts;
    vector<vector<float> > x;
    vector<vector<float> > y;
    vector<vector<float> > z;

    H5ParticleReader::readHDF5("/tmp/PT/points.h5", &ts, &x, &y, &z);

    printf("%d %d %d %d %d %d\n", ts.size(), ts.front().size(), x.size(), x.front().size(), y.size(), z.size());

    return 0;
}
