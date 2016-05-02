
#include "H5ParticleReader.hpp"

#define FILE "/tmp/PT/points.h5"
#define MAX_NAME 256

using namespace std;

void H5ParticleReader::extract1DFloatData(hid_t group, char * datasetName, vector<float> * vData )
{
    int i, ndims;
    hid_t dsid, sid, tid;
    herr_t status;

    hsize_t * ndimsT;

    hsize_t size;

    if( vData == NULL )
    {
        return;
    }

    dsid = H5Dopen(group, datasetName, H5P_DEFAULT);

    sid = H5Dget_space(dsid);
    tid = H5Dget_type(dsid);

    ndims = H5Sget_simple_extent_ndims(sid);
    if(ndims != 1)
    { return; }

    ndimsT = new hsize_t[ndims];

    H5Sget_simple_extent_dims(sid, ndimsT, NULL );

    size = ndimsT[0];

    // reserve some space in the vector
    vData->resize(size);

    status = H5Dread(dsid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vData->data());

    delete[] ndimsT;

    H5Sclose(sid);
    H5Tclose(tid);

    status = H5Dclose(dsid);
}

/* Dipatches data initially packed in vX as x, y, z, x, y, z ...
 * By moving y coordinates to vY vector and z coordinates in vZ vector
 */
void H5ParticleReader::dispatch(vector<float> * vX, vector<float> * vY, vector<float> * vZ)
{
    int idx = 1;
    while(idx < vX->size())
    {
        vY->push_back(vX->at(idx));
        vX->erase(vX->begin() + idx);
        idx = idx + 2;
    } 
    idx = 1;
    while(idx < vX->size())
    {
        vZ->push_back(vX->at(idx));
        vX->erase(vX->begin() + idx);
        idx = idx + 1;
    } 
}

void H5ParticleReader::fillTimeSteps(vector<float> * vX, vector<float> * vTimeStep)
{
    for(int i = 1; i < vX->size(); i++)
    {
        vTimeStep->push_back(vTimeStep->back() + 1.0);
    } 
}

int H5ParticleReader::readHDF5(char * filename, vector<vector<float> > * vvTimeStep, vector< vector<float> > * vvX, vector< vector<float> > * vvY, vector< vector<float> > * vvZ)
{
    int i, j, k;
    herr_t      status;
    char name[MAX_NAME];

    /* Ensure that all the pointers are correctly passed */
    if(!vvTimeStep || !vvX || !vvY || !vvZ)
    {
       return 1; 
    }

    /* Open an existing file. */
    hid_t fileID = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t rootGID = H5Gopen(fileID, "/", H5P_DEFAULT);

    /*status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, */
            /*dset_data);*/

    hsize_t nbPoints;
    status = H5Gget_num_objs(rootGID, &nbPoints);
    std::cout << "Number of Points: " << nbPoints << std::endl;

    for(i = 0 ; i < nbPoints; i++)
    {
        vector<float> vTimeStep;
        vector<float> vX; 
        vector<float> vY; 
        vector<float> vZ;

        hssize_t len = H5Gget_objname_by_idx(rootGID, (hsize_t)i, name, (size_t)MAX_NAME );

        hid_t currentPointGID = H5Gopen(rootGID, name, H5P_DEFAULT);

        hsize_t nbObjects;
        status = H5Gget_num_objs(currentPointGID, &nbObjects);
        //std::cout << "Opening group " << name << " (nbObjects=" << nbObjects << ")" << std::endl;

        for(j = 0 ; j < nbObjects; j++)
        {
            hssize_t len = H5Gget_objname_by_idx(currentPointGID, (hsize_t)j, name, (size_t)MAX_NAME );
            //std::cout << "Object " << name << std::endl;

            if(strcmp(name, "starttime") == 0)
            {
                extract1DFloatData(currentPointGID, name, &vTimeStep);
                float * data = (float *)(vTimeStep.data());

                std::cout << "Start time: " << data[0] << std::endl;
            }
            if(strcmp(name, "coords") == 0)
            {
                extract1DFloatData(currentPointGID, name, &vX);
                float * data = (float *)(vX.data());
                int size = vX.size();

                std::cout << "Read " << size << " points (Integrity: " << ((size % 3) == 0) << std::endl;
                std::cout << "Extract: ";
                for(k = 0; k < 4; k++)
                {
                    std::cout << "(" << data[k * 3 + 0] << ", " << data[k * 3 + 1] << ", " << data[k * 3 + 2] << ") ; ";
                }
                std::cout << "..." << std::endl;
            }
        }

        /* dispatch (x,y,z);(x,y,z) ... packed in vX in other arrays */
        dispatch(&vX, &vY, &vZ);

        fillTimeSteps(&vX, &vTimeStep);

        float * data = (float *)(vX.data());
        std::cout << "Extract: vX";
        for(k = 0; k < 4; k++)
        {
            std::cout << data[k] << " ";
        }
        std::cout << "..." << std::endl;
        data = (float *)(vY.data());
        std::cout << "Extract: vY";
        for(k = 0; k < 4; k++)
        {
            std::cout << "(" << data[k * 3 + 0] << ", " << data[k * 3 + 1] << ", " << data[k * 3 + 2] << ") ; ";
        }
        std::cout << "..." << std::endl;
        data = (float *)(vZ.data());
        std::cout << "Extract: vZ";
        for(k = 0; k < 4; k++)
        {
            std::cout << "(" << data[k * 3 + 0] << ", " << data[k * 3 + 1] << ", " << data[k * 3 + 2] << ") ; ";
        }
        std::cout << "..." << std::endl;

        /* Push back the built vector */
        std::cout << "Adding vector" << std::endl;
        vvTimeStep->push_back(vTimeStep);
        vvX->push_back(vX);
        vvY->push_back(vY);
        vvZ->push_back(vZ);
        std::cout << "Post-Adding vector" << std::endl;

        status = H5Gclose(currentPointGID);
    }

    /* Close the dataset. */
    //status = H5Dclose(dataset_id);

    status = H5Gclose(rootGID);

    /* Close the file. */
    status = H5Fclose(fileID);

    return 0;
}

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
