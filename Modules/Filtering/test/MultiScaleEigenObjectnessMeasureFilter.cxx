/*=========================================================================
 *
 *  Copyright VivaBrain Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <string>
#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <list>
#include "itkCastImageFilter.h"

#include "itkMultiScaleEigenObjectnessFilter3D.h"
#include "itkMultiScaleEigenObjectnessFilter2D.h"
#include "itkEigenToObjectnessFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include <itkMinimumMaximumImageCalculator.h>
#include "itkBinaryThresholdImageFilter.h"
#include <itksys/SystemTools.hxx>

#define EPSILON  1e-03

int EigenToObjectnessMeasure2D(std::string imageFilename, float sMin, float sMax, int nSteps, int gamma, float tMin, float tMax, bool ScaleObj)
{
  bool objBright = true;

  // Define the dimension of the images
  const unsigned char Dim = 2;

  typedef float             InputPixelType;
  typedef float             OutputPixelType;
  typedef unsigned char     UcharPixelType;
        
  typedef itk::Vector<InputPixelType, Dim>       VectorPixelType;
  // Declare the types of the images
  typedef itk::Image<InputPixelType,Dim>         InputImageType;
  typedef itk::Image<OutputPixelType,Dim>        OutputImageType;
  typedef itk::Image<UcharPixelType,Dim>         UcharImageType;
  typedef itk::Image<VectorPixelType, Dim>       VectorImageType;

  typedef itk::ImageFileReader<InputImageType>   FileReaderType;
  typedef itk::ImageFileWriter<OutputImageType>  FileWriterType;
  typedef itk::ImageFileWriter<UcharImageType>   UcharWriterType;

  typedef itk::RescaleIntensityImageFilter<OutputImageType, UcharImageType> RescaleFilterType;      

  // Declare the type of enhancement filter
  typedef itk:: EigenToObjectnessFilter <InputPixelType, OutputPixelType, Dim>   ObjectnessFilterType;

  //  Declare the type of multiscale enhancement filter
  typedef itk::MultiScaleEigenObjectnessFilter2D <InputImageType, ObjectnessFilterType, OutputImageType >
        MultiScaleEnhancementFilterType;

  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName(imageFilename.c_str());
  try
    { 
    imageReader->Update();
    }
  catch (itk::ExceptionObject &ex)
    { 
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }
        
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput(imageReader->GetOutput());

  ObjectnessFilterType* objectnessFilter = multiScaleEnhancementFilter->GetHessianToMeasureFilter();
  objectnessFilter->SetScaleObjectnessMeasure(ScaleObj);
  objectnessFilter->SetBrightObject(objBright);
  objectnessFilter->SetGamma(gamma);
  objectnessFilter->SetObjectDimension( 1 );

  multiScaleEnhancementFilter->SetSigmaMin(sMin);
  multiScaleEnhancementFilter->SetSigmaMax(sMax);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nSteps);
  
  multiScaleEnhancementFilter->SetUpperThreshold( tMax );
  multiScaleEnhancementFilter->SetLowerThreshold( tMin );

  try
    {
    multiScaleEnhancementFilter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  OutputImageType::Pointer  objResult = multiScaleEnhancementFilter->GetVesselnessOutput();       
  OutputImageType::Pointer  vec1Result = multiScaleEnhancementFilter->GetEigenVector1Output();
  OutputImageType::Pointer  vec2Result = multiScaleEnhancementFilter->GetEigenVector2Output();
  
  OutputImageType::Pointer  RBResult = multiScaleEnhancementFilter->GetRBOutput();
  OutputImageType::Pointer  RCResult = multiScaleEnhancementFilter->GetRCOutput();

  OutputImageType::Pointer  scaleResult = multiScaleEnhancementFilter->GetScalesOutput();

  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput(objResult);
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);

  UcharWriterType::Pointer ewriter = UcharWriterType::New();
  

  std::string fwe = itksys::SystemTools::GetFilenameWithoutExtension(imageFilename);
  
  std::stringstream filename;
  filename << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << ".mha";

  ewriter->SetFileName(filename.str().c_str());
  ewriter->SetInput(rescale->GetOutput());
  ewriter->SetUseCompression(true);
  try
    {
    ewriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing vesselness "<< e << std::endl;
    }
        
  if(tMin && tMax)
    {
    OutputImageType::Pointer  thresholdResult = multiScaleEnhancementFilter->GetHystheresisThresholdOutput();
       
    FileWriterType::Pointer twriter = FileWriterType::New();
    
    std::stringstream tfilename;
    tfilename << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_T_" << tMin << "-" << tMax << ".mha";
    twriter->SetFileName(tfilename.str().c_str());
    twriter->SetInput(thresholdResult);
    twriter->SetUseCompression(true);
    try
      {
      twriter->Update();
      }
    catch (itk::ExceptionObject &e)
      {
      std::cerr << e << std::endl;
      }
    }

  FileWriterType::Pointer ewriter1 = FileWriterType::New();
  std::stringstream filename1;
  filename1 << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_E1.mha";
  ewriter1->SetFileName(filename1.str().c_str());
  ewriter1->SetInput(vec1Result);
  try
    {
    ewriter1->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << "Writing vec1 "<<e << std::endl;
    }

  FileWriterType::Pointer ewriter2 = FileWriterType::New();
  std::stringstream filename2;
  filename2 << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_E2.mha";
  ewriter2->SetFileName(filename1.str().c_str());
  ewriter2->SetInput(vec2Result);
  ewriter2->SetUseCompression(true);
  try
    {
    ewriter2->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << "Writing vec2 "<<e << std::endl;
    }

  // Write the image containing the best response scales
  FileWriterType::Pointer writerScales = FileWriterType::New();
  std::stringstream filenamescales;
  filenamescales << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_Scales.mha";
  writerScales->SetFileName(filenamescales.str().c_str());
  writerScales->SetInput(scaleResult);
  try
    {
    writerScales->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << "Writing scales "<<e << std::endl;
    }

  FileWriterType::Pointer writerRB = FileWriterType::New();
  std::stringstream filenameRB;
  filenameRB << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_RB.mha";
  writerRB->SetFileName(filenameRB.str().c_str());
  writerRB->SetInput(RBResult);
  writerRB->SetUseCompression(true);
  try
    {
    writerRB->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing RB "<< e << std::endl;
    }
  FileWriterType::Pointer writerRC = FileWriterType::New();
  std::stringstream filenameRC;
  filenameRC << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_RC.mha";
  writerRC->SetFileName(filenameRC.str().c_str());
  writerRC->SetInput(RCResult);
  writerRC->SetUseCompression(true);
  try
    {
    writerRC->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing RC "<< e << std::endl;
    }

  return 0;
};

int EigenToObjectnessMeasure3D(std::string imageFilename, float sMin, float sMax, int nSteps, int gamma, float tMin, float tMax, bool ScaleObj)
{
  // Define the dimension of the images
  const unsigned char Dim = 3;
  bool objBright = true;

  typedef float             InputPixelType;
  typedef float             OutputPixelType;
  typedef unsigned char     UcharPixelType;

  typedef itk::Vector<InputPixelType, Dim>    VectorPixelType;
  typedef itk::Image<InputPixelType,Dim>      InputImageType;  
  typedef itk::Image<OutputPixelType,Dim>     OutputImageType;
  typedef itk::Image<UcharPixelType,Dim>      UcharImageType;
  typedef itk::Image<VectorPixelType, Dim>    VectorImageType;

  typedef itk::ImageFileReader <InputImageType>   FileReaderType;
  typedef itk::ImageFileWriter<OutputImageType>   FileWriterType;
  typedef itk::ImageFileWriter<UcharImageType>    UcharWriterType;
  
  typedef itk::RescaleIntensityImageFilter<OutputImageType, UcharImageType> RescaleFilterType;
        
  // Declare the type of enhancement filter        
  typedef itk::EigenToObjectnessFilter <InputPixelType, OutputPixelType, Dim>   ObjectnessFilterType;

  //  Declare the type of multiscale enhancement filter
  typedef itk::MultiScaleEigenObjectnessFilter3D< InputImageType, ObjectnessFilterType, OutputImageType >  MultiScaleEnhancementFilterType;

  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName(imageFilename.c_str());
  try
    { 
    imageReader->Update();
    }
  catch (itk::ExceptionObject &ex)
    { 
    std::cout << ex << std::endl;
    }

  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput(imageReader->GetOutput());

  ObjectnessFilterType* objectnessFilter = multiScaleEnhancementFilter->GetHessianToMeasureFilter();
  objectnessFilter->SetScaleObjectnessMeasure(ScaleObj);
  objectnessFilter->SetBrightObject(objBright);
  objectnessFilter->SetAlpha(0.25);
  objectnessFilter->SetBeta(0.25);
  objectnessFilter->SetGamma(gamma);
  objectnessFilter->SetObjectDimension( 1 );

  multiScaleEnhancementFilter->SetSigmaMin(sMin);
  multiScaleEnhancementFilter->SetSigmaMax(sMax);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nSteps);
  multiScaleEnhancementFilter->SetUpperThreshold(tMax);
  multiScaleEnhancementFilter->SetLowerThreshold(tMin);  

  try
    {
    multiScaleEnhancementFilter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  OutputImageType::Pointer  objResult = multiScaleEnhancementFilter->GetVesselnessOutput();

  OutputImageType::Pointer  thresholdResult = multiScaleEnhancementFilter->GetHystheresisThresholdOutput();
  
  OutputImageType::Pointer  vec1Result = multiScaleEnhancementFilter->GetEigenVector1Output();
  OutputImageType::Pointer  vec2Result = multiScaleEnhancementFilter->GetEigenVector2Output();
  OutputImageType::Pointer  vec3Result = multiScaleEnhancementFilter->GetEigenVector3Output();

  OutputImageType::Pointer  RAResult = multiScaleEnhancementFilter->GetRAOutput();
  OutputImageType::Pointer  RBResult = multiScaleEnhancementFilter->GetRBOutput();
  OutputImageType::Pointer  RCResult = multiScaleEnhancementFilter->GetRCOutput();

  OutputImageType::Pointer  scaleResult = multiScaleEnhancementFilter->GetScalesOutput();
  
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput(objResult);
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);

  UcharWriterType::Pointer ewriter = UcharWriterType::New();

  std::string fwe = itksys::SystemTools::GetFilenameWithoutExtension(imageFilename);
  std::string path = itksys::SystemTools::GetFilenamePath(imageFilename);
  
  std::stringstream filename;
  filename << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << ".mha";

  ewriter->SetFileName(filename.str().c_str());
  ewriter->SetInput(rescale->GetOutput());
  ewriter->SetUseCompression(true);
  try
    {
    ewriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing vesselness "<< e << std::endl;
    }
        
  if(tMin && tMax)
    {
    FileWriterType::Pointer twriter = FileWriterType::New();
    std::stringstream tfilename;
    tfilename << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_T_" << tMin << "-" << tMax << ".mha";
    twriter->SetFileName(tfilename.str().c_str());
    twriter->SetInput(thresholdResult);
    try
      {
      twriter->Update();
      }
    catch (itk::ExceptionObject &e)
      {
      std::cerr << e << std::endl;
      }
    }

  FileWriterType::Pointer ewriter1 = FileWriterType::New();
  std::stringstream filename1;
  filename1 << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_E1.mha";
  ewriter1->SetFileName(filename1.str().c_str());
  ewriter1->SetInput(vec1Result);
  ewriter1->SetUseCompression(true);
  try
    {
    ewriter1->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer ewriter2 = FileWriterType::New();
  std::stringstream filename2;
  filename2 << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_E2.mha";
  ewriter2->SetFileName(filename2.str().c_str());
  ewriter2->SetInput(vec2Result);
  ewriter2->SetUseCompression(true);
  try
    {
    ewriter2->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer ewriter3 = FileWriterType::New();
  std::stringstream filename3;
  filename3 << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_E3.mha";
  ewriter3->SetFileName(filename3.str().c_str());
  ewriter3->SetInput(vec3Result);
  ewriter3->SetUseCompression(true);
  try
    {
    ewriter3->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  // Write the image containing the best response scales
  FileWriterType::Pointer writerScales = FileWriterType::New();
  std::stringstream filenameScales;
  filenameScales << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_Scales.mha";
  writerScales->SetFileName(filenameScales.str().c_str());
  writerScales->SetInput(scaleResult);
  writerScales->SetUseCompression(true);
  try
    {
    writerScales->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer writerRA = FileWriterType::New();
  std::stringstream filenameRA;
  filenameRA << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_RA.mha";
  writerRA->SetFileName(filenameRA.str().c_str());
  writerRA->SetInput(RAResult);
  writerRA->SetUseCompression(true);
  try
    {
    writerRA->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer writerRB = FileWriterType::New();
  std::stringstream filenameRB;
  filenameRB << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_RB.mha";
  writerRB->SetFileName(filenameRB.str().c_str());
  writerRB->SetInput(RBResult);
  writerRB->SetUseCompression(true);
  try
    {
    writerRB->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }
  FileWriterType::Pointer writerRC = FileWriterType::New();
  std::stringstream filenameRC;
  filenameRC << fwe << "_" << sMin << "-" << sMax << "_" << nSteps << "_RC.mha";
  writerRC->SetFileName(filenameRC.str().c_str());
  writerRC->SetInput(RCResult);
  writerRC->SetUseCompression(true);
  try
    {
    writerRC->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  typedef itk::BinaryThresholdImageFilter< OutputImageType, UcharImageType >  FilterType;
  FilterType::Pointer filter = FilterType::New();
  UcharWriterType::Pointer thresh = UcharWriterType::New();
        
  filter-> SetInput( multiScaleEnhancementFilter->GetScalesOutput());
  double sNew = 0;
  double s = sMin;
        
  for (s = sMin; s <= sMax; s++)
    {
    filter->SetOutsideValue( 0 );
    filter->SetInsideValue( 1 );            
    filter->SetLowerThreshold( s );
    sNew = s + 1;
    filter->SetUpperThreshold( sNew );
    filter->Update();
            
    thresh->SetInput(filter->GetOutput());
    
    std::cout << "s "<< s << " sNew " << sNew << " sMin " << sMin << " sMax " << sMax << std::endl;
     
    std::stringstream filenameThresh;
    filenameThresh << fwe << "_s_" << s << "-" << sNew << ".mha";
    thresh->SetFileName(filenameThresh.str().c_str());
    thresh->SetUseCompression(true);
    thresh->Update();
    }
  return 0;
};

/** Main test function */
int MultiScaleEigenObjectnessMeasureFilter(int argc, char* argv[])
{
  if(argc < 9)
    {
    std::cerr << "Missing parameters: " << argv[0] << " image dimension(2/3) sMin sMax nSteps gamma [Tmin Tmax ScaleObjectness]" << std::endl;
    return EXIT_FAILURE;
    }  
  
  std::string img = argv[1]; 
  int dimension = atoi(argv[2]);
  float sMin = atof(argv[3]);
  float sMax = atof(argv[4]);
  int nSteps = atoi(argv[5]);
  int gamma = atoi(argv[6]);
  float Tmin, Tmax;     
  bool SO=false;
        
  if(argv[7] && argv[8])
    {
    Tmin = atof(argv[7]);
    Tmax = atof(argv[8]);
    }
  else
    {
    Tmin = 0;
    Tmax = 0;
    }
        
  if(strcmp(argv[9],"1"))
    {
    SO=true;
    }
  else
    {
    SO=false;
    }        
  int r = 0;  

  if(dimension == 2)
    {
    r = EigenToObjectnessMeasure2D(img, sMin, sMax, nSteps, gamma, Tmin, Tmax, SO);
    }
  else if(dimension == 3)
    {
    r = EigenToObjectnessMeasure3D(img, sMin, sMax, nSteps, gamma, Tmin, Tmax, SO);
    }
  else
    {
    return EXIT_FAILURE;
    }
  return 0;
}
