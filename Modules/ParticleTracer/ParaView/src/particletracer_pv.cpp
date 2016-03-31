
#include <vtkDataArrayCollection.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkEnSightGoldReader.h>
#include <vtkInformation.h>
#include <vtkParticleTracer.h>
#include <vtkPointSource.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkPVOptions.h>
#include <vtkInitializationHelper.h>
#include <vtkProcessModule.h>

#include <vtksys/SystemTools.hxx>


#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

using namespace boost;
namespace po = boost::program_options;

#ifdef _cplusplus
extern "C" {
#endif
#include "hdf5.h"
#ifdef _cplusplus
}
#endif

hid_t fileID;
hid_t plist;

/* Initialize hdf5 reading */
void initHDF5(const char * filename)
{
    /* Open the file (Truncate it if it already exists */
    fileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* To enable extendable datasets, we need to enable chunking */
    /* So do the following lines */
    plist = H5Pcreate (H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] = {500};
    H5Pset_layout( plist, H5D_CHUNKED ); 
    herr_t status = H5Pset_chunk (plist, 1, chunk_dims);
}

/* Finalize hdf5 reading */
void finalizeHDF5()
{
    H5Fclose(fileID);
}

int main(int argc, char ** argv)
{
    /* Parses arguments using boost::program_options */
    po::variables_map vm;

    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("input.data.filename", po::value<std::string>(), "Sets the hdf5 output filename.")
            ("input.pointsources.filename", po::value<std::string>()->default_value("points.csv"), "Sets the input point source csv filename (format: pos_x,pos_y,pos_z,radius,nbpoints).")
            ("output.filename", po::value<std::string>()->default_value("/tmp/PT/points.h5"), "Sets the hdf5 output filename.")
            ("output.vtk.enable", po::value<bool>()->default_value(false), "Sets whether to output vtk files or not.")
            ;

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    /* Read the info about the point sources */
    /* Data is read from a csv file */
    /* Each line represents: */
    /* position_x, position_y, position_z, radius, number_of_points */
    std::vector<vtkSmartPointer<vtkPointSource> > vps;
    if(vm.count("input.pointsources.filename"))
    {
        /* Open the file given as argument */
        std::ifstream infile(vm["input.pointsources.filename"].as<std::string>().c_str());
        std::string line;
        /* Iterate over the lines in the file */
        while (std::getline(infile, line))
        {
            /* Tokenize the line and store the result in result */
            std::vector<std::string> result;
            typedef boost::tokenizer< boost::escaped_list_separator<char> , std::string::const_iterator, std::string> Tokenizer;
            Tokenizer tok(line);
            result.assign(tok.begin(),tok.end());

            /* Check that we read enough elements */
            if(result.size() < 5)
            {
                std::cout << "Bad entry: " << line << std::endl;
                return 0;
            }

            /* Create a point source for the particle tracer */
            vtkSmartPointer<vtkPointSource> ps = vtkSmartPointer<vtkPointSource>::New();
            ps->SetCenter(boost::lexical_cast<double>(result[0]), boost::lexical_cast<double>(result[1]), boost::lexical_cast<double>(result[2]));
            ps->SetRadius(boost::lexical_cast<double>(result[3]));
            ps->SetNumberOfPoints(boost::lexical_cast<int>(result[4]));

            /* Store it in a vector */
            vps.push_back(ps);
        }
        /* close the file */
        infile.close();
    }
    else
    {
        std::cout << "Input points file was not set (--input.pointsources.filename)" << std::endl;
        return 1;
    }

    /* If we don't have any vtkPointSource data, we exit the application */
    if(vps.size() == 0)
    {
        std::cout << "No input point sources were set (--input.pointsources.filename)" << std::endl;
        return 1;
    }

    /* Initialize ParaView environement */
    /* Needed for the particle tracer */
    static vtkSmartPointer<vtkPVOptions> options = vtkSmartPointer<vtkPVOptions>::New();
    char** execName = new char*[1];
    execName[0] = vtksys::SystemTools::DuplicateString("foobar_string");
    vtkInitializationHelper::Initialize( 1, execName, vtkProcessModule::PROCESS_CLIENT, options );

    /* Read input dataset (Ensight gold format for now) */
    if(! vm.count("input.data.filename"))
    {
        std::cout << "No input data specified (--input.data.filename)" << std::endl;
        return 1;
    }
    vtkSmartPointer<vtkEnSightGoldBinaryReader> egr = vtkEnSightGoldBinaryReader::New();
    //egr->SetCaseFileName("/ssd/ancel/angiotk/pipelines/script_angiotk_phantom/feelpp/export/applications/models/fluid/phantom/np_16/fluid.exports/Export.case");
    egr->SetCaseFileName(vm["input.data.filename"].as<std::string>().c_str());
    egr->Update();

    /* Get Information about input data */
    std::cout << egr->GetTimeValue() << " " << egr->GetMaximumTimeValue() << std::endl;

    vtkDataArrayCollection* dac = egr->GetTimeSets();
    std::cout << dac->GetNumberOfItems() << std::endl;
    dac->InitTraversal();
    vtkDataArray * da = dac->GetNextItem();
    std::cout << da->GetNumberOfTuples() << std::endl;
    //for(int i = 0; i < da->GetNumberOfTuples(); i=i+10)
    //{
        //std::cout << da->GetTuple1(i) << std::endl;
    //}

    std::cout << "Using " << vps.size() << " point sources" << std::endl;

    /* Initialize the particle tracer */
    vtkSmartPointer<vtkParticleTracer> pt = vtkSmartPointer<vtkParticleTracer>::New();
    pt->SetInputConnection(0, egr->GetOutputPort());
    for(int i = 0; i < 1; i++)
    {
        pt->SetInputConnection(i + 1, vps[i]->GetOutputPort());
    }
    pt->SetComputeVorticity(0);
    pt->SetForceReinjectionEveryNSteps(12);

    /* Initialize output file */
    if(vm.count("output.filename"))
    { initHDF5(vm["output.filename"].as<std::string>().c_str()); }
    else
    {
        std::cout << "Output hdf5 file was not set (--output.filename)" << std::endl;
        return 1;
    }

    hsize_t entryDims[1] = {3};
    hsize_t offset[1] = {0};
    hsize_t maxDims[1] = {H5S_UNLIMITED};

    hid_t rootGID = H5Gopen(fileID, "/", H5P_DEFAULT);

    /* Iterate over all the timesteps */
    double t1 = da->GetTuple1(0); 
    double t2;
    for(int i = 1; i < da->GetNumberOfTuples(); i++)
    {
        /* Display time information */
        t2 = da->GetTuple1(i);
        t1 = da->GetTuple1(0);
        std::cout << t1 << " -> " << t2 << std::endl;

        pt->SetStartTime(t1);
        pt->SetTerminationTime(t2);
        pt->Update();

        /* Get particle data */
        vtkPolyData * pd = pt->GetOutput();
        //std::cout << pd->GetPoints()->GetNumberOfPoints() << std::endl;

        /* We can also access useful field data in the data */
        vtkDataArray * pids = pd->GetPointData()->GetScalars("ParticleId");
        //if(pids)
        //{
            //std::cout << pids->GetNumberOfComponents() << " " << pids->GetNumberOfTuples() << std::endl;
        //}
        vtkDataArray * page = pd->GetPointData()->GetScalars("ParticleAge");
        
        /* Enable outputting data to vtk format if required */
        /* Allow the users to check that the particle tracer behaves accordingly */
        if(vm.count("output.vtk.enable"))
        {
            if(vm["output.vtk.enable"].as<bool>())
            {
                std::ostringstream oss;
                oss << "particles." << i << ".vtp";
                vtkSmartPointer<vtkXMLPolyDataWriter> pdw = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
                    pdw->SetInputData(pt->GetOutput()); 
                    pdw->SetFileName(oss.str().c_str());
                    pdw->Write();
            }
        }

        /* Iterate over the points (each particle) and save the associated data */
        for(int j = 0; j < pd->GetNumberOfPoints(); j++)
        {
            std::ostringstream oss;
            oss << j << std::endl;

            /* Test if the group corresponding to the particle id under the root node exists */
            /* Create it if not */
            hid_t groupID;
            if( H5Lexists( rootGID, oss.str().c_str(), H5P_DEFAULT ) )
            { groupID = H5Gopen(rootGID, oss.str().c_str(), H5P_DEFAULT); }
            else
            { groupID = H5Gcreate(rootGID, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); }

            /* Test if the coords dataset exists for the current particle group */
            /* If not, we create it */
            if(!(H5Lexists( groupID, "coords", H5P_DEFAULT )))
            {
                /* Create a dataspace to store data */
                /* Only store a 3-uple corresponding to the first (x, y, z) coordinates */
                hid_t dataspace = H5Screate_simple(1, entryDims, maxDims); 
                hid_t dataset = H5Dcreate(groupID, "coords", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

                /* Get the current point position and store it in the file */
                double * data = pd->GetPoint(j);
                /* Convert data to float (needed by JEMRIS) */
                float fdata[3];
                for(int i = 0; i < 3; i++)
                { fdata[i] = (float)(data[i]); }
                herr_t status = H5Dwrite (dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fdata);

                /* close all opened objects */
                H5Dclose(dataset);
                H5Sclose(dataspace);
            }
            /* Here we switch to an append mode for the hdf5 file */
            else
            {
                herr_t status;
                /* Open the coords dataset */
                hid_t dataset = H5Dopen(groupID, "coords", H5P_DEFAULT);

                /* Get information about the current dimension of the dataset */
                hsize_t curDims[1];
                hsize_t curMaxDims[1];
                hid_t filespace = H5Dget_space (dataset);
                H5Sget_simple_extent_dims(filespace, curDims, curMaxDims); 
                H5Sclose(filespace);

                /* Fill the offset (from where we will be writing) and the current dimensions (to extend the dataset) */ 
                offset[0] = curDims[0];
                curDims[0] = curDims[0] + entryDims[0];

                /* We extend the dataset and get the new data size */
                status = H5Dset_extent (dataset, curDims);
                filespace = H5Dget_space (dataset);
                
                /* We select the hyperslab corresponding to the current 3-uple that we want to fill */
                /* with teh previoulsy computed offset */
                status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, entryDims, NULL);  
                //std::cout << "offset[0]:" << offset[0] << " " << " curDims[0]: " << curDims[0] << " entryDims[0]: " << entryDims[0] << std::endl;

                hid_t memspace = H5Screate_simple(1, entryDims, NULL); 

                double * data = pd->GetPoint(j);

                /* Convert data to float (needed by JEMRIS) */
                float fdata[3];
                for(int i = 0; i < 3; i++)
                { fdata[i] = (float)(data[i]); }

                /* Write data to the file */
                status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, fdata);

                /* Close all unused object */
                H5Sclose(filespace);
                H5Sclose(memspace);

                H5Dclose(dataset);
            }

            /* Check if the starttime dataset exists */
            if(!(H5Lexists( groupID, "starttime", H5P_DEFAULT )))
            {
                /* Store the start time for the current particle */
                /* (Here we are at the point when we first encounter it */
                hsize_t st_dims[1] = {1};
                hsize_t st_max_dims[1] = {1};
                hid_t dataspace = H5Screate_simple(1, st_dims, st_max_dims); 
                hid_t dataset = H5Dcreate(groupID, "starttime", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                int starttime = 0;
                herr_t status = H5Dwrite (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &starttime);
                H5Dclose(dataset);
            }
            
            /* Close the current particle group */
            H5Gclose(groupID);
        }

        t1 = t2;
    }

    finalizeHDF5();

#if 0
    vtkSmartPointer<vtkPolyDataMapper> pdm = vtkSmartPointer<vtkPolyDataMapper>::New();
        pdm->SetInputConnection(pt->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(pdm);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    renderer->AddActor(actor);
      renderer->SetBackground(1,1,1); //white
       
        // Render and interact
           renderWindow->Render();
             renderWindowInteractor->Start();
#endif

    //for(int i = 0; i < da->GetNumberOfTuples(); i=i+10)
    //{
        //vtkInformation * info = pt->GetInformation();

        //info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), da->GetTuple1(i));

        //pt->Update();
    //}

	return 0;
}
