// Need to add some copyright header

//#include <valarray>
#include <string>
//#include <vector>
#include <iostream>

//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <list>
#include "itkCastImageFilter.h"
//#include "vnl/vnl_math.h"

#include "itkMultiScaleEigenObjectnessFilter3D.h"
#include "itkMultiScaleEigenObjectnessFilter2D.h"
#include "itkEigenToObjectnessFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include <itkMinimumMaximumImageCalculator.h>
#include "itkBinaryThresholdImageFilter.h"

#define EPSILON  1e-03

int EigenToObjectnessMeasure2D(std::string file, std::string fe, std::string path, float smin, float smax, int ns, int gamma, float Tmin, float Tmax, bool SO)
{
  std::string  img = file;
  std::string   folder = path;
  std::string ext = fe;

  float sMin = smin;
  float sMax = smax;
  int nSteps = ns; 
        
  bool objBright = true;
  bool ScaleObj = SO;

  std::stringstream osMin;
  std::stringstream osMax;
  std::stringstream onSteps;
  std::stringstream tMin;
  std::stringstream tMax;
  osMin<<sMin;
  osMax<<sMax;
  onSteps<<nSteps;
        
  tMin<<Tmin;
  tMax<<Tmax;

  // Define the dimension of the images
  const unsigned char Dim = 2;

  typedef float             InputPixelType;
  typedef float            OutputPixelType;
  typedef unsigned char          UcharPixelType;
        
  typedef itk::Vector<InputPixelType, Dim>      VectorPixelType;
  // Declare the types of the images
  typedef itk::Image<InputPixelType,Dim>        InputImageType;
  typedef itk::Image<OutputPixelType,Dim>        OutputImageType;
  typedef itk::Image<UcharPixelType,Dim>        UcharImageType;
  typedef itk::Image<VectorPixelType, Dim>      VectorImageType;

  typedef itk::ImageFileReader<InputImageType>    FileReaderType;
  typedef itk::ImageFileWriter<OutputImageType>    FileWriterType;
  typedef itk::ImageFileWriter<UcharImageType>    UcharWriterType;
  typedef itk::ImageFileWriter<VectorImageType>    VectorWriterType;

  typedef itk::RescaleIntensityImageFilter<OutputImageType, UcharImageType> RescaleFilterType;      
  typedef itk::CastImageFilter< OutputImageType, UcharImageType> CastToRealFilterType;

  // Declare the type of enhancement filter
  typedef itk:: EigenToObjectnessFilter <InputPixelType, OutputPixelType, Dim>   ObjectnessFilterType;

  //  Declare the type of multiscale enhancement filter
  typedef itk::MultiScaleEigenObjectnessFilter2D <InputImageType, ObjectnessFilterType, OutputImageType >
        MultiScaleEnhancementFilterType;

  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName(folder + img+ ext);
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
  multiScaleEnhancementFilter->SetSigmaMax( sMax);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nSteps);
  
  multiScaleEnhancementFilter->SetUpperThreshold( Tmax );
  multiScaleEnhancementFilter->SetLowerThreshold( Tmin );

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
  ewriter->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + ".mhd");
  ewriter->SetInput(rescale->GetOutput());
  try
    {
    ewriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing vesselness "<< e << std::endl;
    }
        
  if(Tmin && Tmax)
    {
    OutputImageType::Pointer  thresholdResult = multiScaleEnhancementFilter->GetHystheresisThresholdOutput();
       
    FileWriterType::Pointer twriter = FileWriterType::New();
    twriter->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_T_" + tMin.str() + "-" + tMax.str() + ".mhd");
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
  ewriter1->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_E1" + ".mhd");
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
  ewriter2->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_"+ onSteps.str() + "_E2" + ".mhd");
  ewriter2->SetInput(vec2Result);
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
  writerScales->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_"+ onSteps.str() + "_Scales" + ".mhd");
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
  writerRB->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_RB"+ ".mhd");
  writerRB->SetInput(RBResult);
  try
    {
    writerRB->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing RB "<< e << std::endl;
    }
  FileWriterType::Pointer writerRC = FileWriterType::New();
  writerRC->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_RC" + ".mhd");
  writerRC->SetInput(RCResult);
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

int EigenToObjectnessMeasure3D(std::string file, std::string fe, std::string path, float smin, float smax, int ns, int gamma, float Tmin, float Tmax, bool SO)
{
  std::string  img = file;
  std::string   folder = path;
  std::string ext = fe;

  float sMin = smin;
  float sMax = smax;
  int nSteps = ns; 
        
  bool objBright = true;
  bool ScaleObj = SO;

  std::stringstream osMin;
  std::stringstream osMax;
  std::stringstream onSteps;
  std::stringstream tMin;
  std::stringstream tMax;
  osMin<<sMin;
  osMax<<sMax;
  onSteps<<nSteps;
  tMin<<Tmin;
  tMax<<Tmax;

  // Define the dimension of the images
  const unsigned char Dim = 3;

  typedef float             InputPixelType;
  typedef float            OutputPixelType;
  typedef unsigned char          UcharPixelType;
  typedef itk::Vector<InputPixelType, Dim>    VectorPixelType;
  typedef itk::Image<InputPixelType,Dim>      InputImageType;  
  typedef itk::Image<OutputPixelType,Dim>     OutputImageType;
  typedef itk::Image<UcharPixelType,Dim>      UcharImageType;
  typedef itk::Image<VectorPixelType, Dim>    VectorImageType;

  typedef itk::ImageFileReader <InputImageType>   FileReaderType;
  typedef itk::ImageFileWriter<OutputImageType>   FileWriterType;
  typedef itk::ImageFileWriter<UcharImageType>    UcharWriterType;
  typedef itk::ImageFileWriter<VectorImageType>   VectorWriterType;
  
  typedef itk::RescaleIntensityImageFilter<OutputImageType, UcharImageType> RescaleFilterType;
        
  typedef itk::CastImageFilter< OutputImageType, UcharImageType> CastToRealFilterType;
        
  // Declare the type of enhancement filter        
  typedef itk::EigenToObjectnessFilter <InputPixelType, OutputPixelType, Dim>   ObjectnessFilterType;

  //  Declare the type of multiscale enhancement filter
  typedef itk::MultiScaleEigenObjectnessFilter3D< InputImageType, ObjectnessFilterType, OutputImageType >  MultiScaleEnhancementFilterType;

  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName(folder + img+ ext);
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
  multiScaleEnhancementFilter->SetSigmaMax( sMax);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nSteps);
  multiScaleEnhancementFilter->SetUpperThreshold( Tmax);
  multiScaleEnhancementFilter->SetLowerThreshold( Tmin );  

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
  ewriter->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + ".mhd");
  ewriter->SetInput(rescale->GetOutput());
  try
    {
    ewriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr <<"Writing vesselness "<< e << std::endl;
    }
        
  if(Tmin && Tmax)
    {
    FileWriterType::Pointer twriter = FileWriterType::New();
    twriter->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_T_" + tMin.str() + "-" + tMax.str() + ".mhd");
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
  ewriter1->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_E1" + ".mhd");
  ewriter1->SetInput(vec1Result);
  try
    {
    ewriter1->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer ewriter2 = FileWriterType::New();
  ewriter2->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_E2" + ".mhd");
  ewriter2->SetInput(vec2Result);
  try
    {
    ewriter2->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer ewriter3 = FileWriterType::New();
  ewriter3->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_E3" + ".mhd");
  ewriter3->SetInput(vec3Result);
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
  writerScales->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_Scales" + ".mhd");
  writerScales->SetInput(scaleResult);
  try
    {
    writerScales->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer writerRA = FileWriterType::New();
  writerRA->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_RA"+ ".mhd");
  writerRA->SetInput(RAResult);
  try
    {
    writerRA->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }

  FileWriterType::Pointer writerRB = FileWriterType::New();
  writerRB->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_RB" + ".mhd");
  writerRB->SetInput(RBResult);
  try
    {
    writerRB->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }
  FileWriterType::Pointer writerRC = FileWriterType::New();
  writerRC->SetFileName(folder + img + "_" + osMin.str() + "-" + osMax.str() + "_" + onSteps.str() + "_RC" + ".mhd");
  writerRC->SetInput(RCResult);
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
            
    std::stringstream os;
    std::stringstream osNew;
    os<<s;
    osNew<<sNew;
    std::cout<<"s "<< s<<" sNew "<<sNew<<" sMin "<<sMin<<" sMax "<<sMax<<std::endl;
    thresh->SetFileName( folder + img + "_s_" + os.str() + "-" + osNew.str() + ".mhd");
    thresh->Update();
    }
  return 0;
};

/** Main test function */
int main(int argc, char* argv[])
{
  if(argc < 9)
    {
    std::cerr<<"Missing parameters: "<<argv[0]<<" folder image ext dimension(2/3) sMin sMax nSteps gamma [Tmin Tmax ScaleObjectness]"<<std::endl;
    return EXIT_FAILURE;
    }  
  
  std::string folder = argv[1];
  std::string  img = argv[2];
  std::string ext = argv[3];  
  int dimension = atoi(argv[4]);
  float sMin = atof(argv[5]);
  float sMax = atof(argv[6]);
  int nSteps = atoi(argv[7]);
  int gamma = atoi(argv[8]);
  float Tmin, Tmax;     
  bool SO=false;
        
  if(argv[9] && argv[10])
    {
    Tmin = atof(argv[9]);
    Tmax = atof(argv[10]);
    }
  else
    {
    Tmin = 0;
    Tmax = 0;
    }
        
  if(strcmp(argv[11],"1"))
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
    r = EigenToObjectnessMeasure2D(img, ext, folder, sMin, sMax, nSteps, gamma, Tmin, Tmax, SO);
    }
  else if(dimension == 3)
    {
    r = EigenToObjectnessMeasure3D(img, ext, folder, sMin, sMax, nSteps, gamma, Tmin, Tmax, SO);
    }
  else
    {
    return EXIT_FAILURE;
    }
  return 0;
}
