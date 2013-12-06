/*=========================================================================

=========================================================================*/
#ifndef __itkEigenToObjectnessFilter_hxx
#define __itkEigenToObjectnessFilter_hxx


#include "itkEigenToObjectnessFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

namespace itk
{

/** Constructor */
template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension >
::EigenToObjectnessFilter()
{
  m_Alpha = 0.25;
  m_Beta = 0.25;
  m_Gamma = 5;

  m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  m_SymmetricEigenValueFilter->SetDimension( ImageDimension );
  m_SymmetricEigenValueFilter->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);

  // instantiate the SymmetricEigenVectorAnalysis filter
  m_EigenVectorAnalysisFilter = EigenVectorAnalysisFilterType::New();
  m_EigenVectorAnalysisFilter->SetDimension( ImageDimension ); 
  m_EigenVectorAnalysisFilter->OrderEigenValuesBy(EigenVectorAnalysisFilterType::FunctorType::OrderByValue);

  m_ScaleObjectnessMeasure = true; 

  // by default extract bright lines (equivalent to vesselness)
  m_ObjectDimension = 1;
  m_BrightObject = true;

  this->ProcessObject::SetNumberOfRequiredOutputs(3);
  this->SetNthOutput( 0, this->MakeOutput( 0 ) );
  this->SetNthOutput( 1, this->MakeOutput( 1 ) );
  this->SetNthOutput( 2, this->MakeOutput( 2 ) );
  this->SetNthOutput( 5, this->MakeOutput( 5 ) );
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
typename EigenToObjectnessFilter<TInputPixel,TOutputPixel,VDimension>::OutputImageType *
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension > 
::GetVesselnessOutput( void )
{
  return dynamic_cast<OutputImageType *>(this->ProcessObject::GetOutput( 0 ) );
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
typename EigenToObjectnessFilter<TInputPixel,TOutputPixel,VDimension>::ArrayOutputImageType *
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension > 
::GetEigenVectorsOutput( void )
{
  return dynamic_cast<ArrayOutputImageType *>(this->ProcessObject::GetOutput( 1 ) );
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
typename EigenToObjectnessFilter<TInputPixel,TOutputPixel,VDimension>::ROutputImageType *
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension > 
::GetROutput( void )
{
  return dynamic_cast<ROutputImageType *>(this->ProcessObject::GetOutput( 2 ) );
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
DataObject::Pointer
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension >
::MakeOutput(unsigned int idx)
{
  DataObject::Pointer output;
  switch (idx)
  {
  case 0:
  output = (OutputImageType::New()).GetPointer();
  break;
  case 1:
  output = (ArrayOutputImageType::New()).GetPointer();
  break;
  case 2:
  output = (ROutputImageType::New()).GetPointer();
  break;
  }
  return output.GetPointer();
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
void
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension > 
::GenerateData()
{
  if (m_ObjectDimension >= ImageDimension)
    {
    throw ExceptionObject(__FILE__, __LINE__,"ObjectDimension must be lower than ImageDimension.",ITK_LOCATION);
    }

  m_SymmetricEigenValueFilter->SetInput( this->GetInput() ); //Eigen value analysis
  m_SymmetricEigenValueFilter->Update();

  typedef typename EigenAnalysisFilterType::OutputImageType EigenValueImageType;
  const typename EigenValueImageType::ConstPointer eigenImage = m_SymmetricEigenValueFilter->GetOutput();

  m_EigenVectorAnalysisFilter->SetInput( this->GetInput() );//Eigen vector analysis
  m_EigenVectorAnalysisFilter->Update();
  const typename AnalysisMatrixImageType::ConstPointer eigenVectorMatrixImage = m_EigenVectorAnalysisFilter->GetOutput();

  // walk the region of eigen values and get the objectness measure
  EigenValueArrayType eigenValues;
  ImageRegionConstIterator<EigenValueImageType> it;
  it = ImageRegionConstIterator<EigenValueImageType>(eigenImage, eigenImage->GetRequestedRegion());

  AnalysisMatrixType eigenMatrix;
  ImageRegionConstIterator< AnalysisMatrixImageType > im;
  im = ImageRegionConstIterator<AnalysisMatrixImageType>(eigenVectorMatrixImage, eigenVectorMatrixImage->GetLargestPossibleRegion()); //Iterator of Hessian vectors

  typename OutputImageType::Pointer objOutput = this->GetVesselnessOutput();
  //typename ArrayOutputImageType::Pointer eigenValuesOutput = this->GetEigenValuesOutput();
  typename ArrayOutputImageType::Pointer eigenVectorsOutput = this->GetEigenVectorsOutput();
  typename ROutputImageType::Pointer ROutput = this->GetROutput();

  this->AllocateOutputs();

  ImageRegionIterator<OutputImageType> objit;
  objit = ImageRegionIterator<OutputImageType>(objOutput, objOutput->GetRequestedRegion());
  ImageRegionIterator<ArrayOutputImageType> vecit;
  vecit = ImageRegionIterator<ArrayOutputImageType>(eigenVectorsOutput, eigenVectorsOutput->GetRequestedRegion());

  ImageRegionIterator<ROutputImageType> Rit;
  Rit = ImageRegionIterator<ROutputImageType>(ROutput, ROutput->GetRequestedRegion());

  objit.GoToBegin();
  vecit.GoToBegin();
  Rit.GoToBegin();

  it.GoToBegin();
  im.GoToBegin();

  FixedArray<EigenValueType, ImageDimension> FinalArray;
  InputIndexType InputIndx;
  RArrayType RArray;
  while (!it.IsAtEnd())
    {
    eigenValues = it.Get(); // Get the eigenvalues, vectors
    eigenMatrix = im.Get();

    InputIndx =it.GetIndex();

    EigenValueArrayType sortedEigenValues = eigenValues;
    AnalysisMatrixType sortedEigenVectors = eigenMatrix;
    bool done = false;
    while (!done)
      {
      done = true;
      for (unsigned int i=0; i<ImageDimension-1; i++) // Sort the eigenvalues by increasing magnitude but retain their sign
        {
        if (vnl_math_abs(sortedEigenValues[i]) > vnl_math_abs(sortedEigenValues[i+1]))
          {
          EigenValueType temp = sortedEigenValues[i+1];
          sortedEigenValues[i+1] = sortedEigenValues[i];
          sortedEigenValues[i] = temp;

          for (unsigned int j=0; j<ImageDimension; j++)
            {
            EigenValueType temp = sortedEigenVectors[j][i+1];
            sortedEigenVectors[j][i+1] = sortedEigenVectors[j][i];
            sortedEigenVectors[j][i] = temp;
            }
          done = false;
          }
        }
      }

    // check whether eigenvalues have the right sign
    bool signConstraintsSatisfied= true;
    for (unsigned int i=m_ObjectDimension; i<ImageDimension; i++)
     
  if (!signConstraintsSatisfied)
    {
    FinalArray.Fill(0.0);
    objit.Set(static_cast< OutputPixelType >(0.0));
    //valit.Set( static_cast< OutputArrayPixelType >(FinalArray));
    vecit.Set( static_cast< OutputArrayPixelType >(FinalArray));
    RArray.Fill(0.0);
    Rit.Set(static_cast< RArrayType >(RArray));

    ++it;
    ++im;

    ++objit;
    //++valit;
    ++vecit;
    ++Rit;
    continue;
    }

  EigenValueArrayType sortedAbsEigenValues;
  for (unsigned int i=0; i<ImageDimension; i++)
    {
    sortedAbsEigenValues[i] = vnl_math_abs(sortedEigenValues[i]);
    }

  // initialize the objectness measure
  double objectnessMeasure = 1.0, rA = 1.0, rB = 1.0;

  // compute objectness from eigenvalue ratios and second-order structureness 
  if (m_ObjectDimension < ImageDimension-1)//for lines in 2D and planes in 2D&3D
    { 
    rA = sortedAbsEigenValues[m_ObjectDimension];
    double rADenominatorBase = 1.0;
    for (unsigned int j = m_ObjectDimension+1; j<ImageDimension; j++)
      {
      rADenominatorBase *= sortedAbsEigenValues[j];
      }
    rA /= pow(rADenominatorBase, 1.0 / (ImageDimension-m_ObjectDimension-1));
    objectnessMeasure *= 1.0 - vcl_exp(- 0.5 * vnl_math_sqr(rA) / vnl_math_sqr(m_Alpha));
    }

  if (m_ObjectDimension > 0)//for 2D&3D blobs, lines in 3D
    {
    rB = sortedAbsEigenValues[m_ObjectDimension-1];
    double rBDenominatorBase = 1.0;
    for (unsigned int j=m_ObjectDimension; j<ImageDimension; j++)
      {
      rBDenominatorBase *= sortedAbsEigenValues[j];
      }
    rB /= pow(rBDenominatorBase, 1.0 / (ImageDimension-m_ObjectDimension));
    objectnessMeasure *= vcl_exp(- 0.5 * vnl_math_sqr(rB) / vnl_math_sqr(m_Beta));
    }
  double frobeniusNorm = 0.0;
  for (unsigned int i=0; i<ImageDimension; i++)
    {
    frobeniusNorm += vnl_math_sqr(sortedAbsEigenValues[i]);
    }
  frobeniusNorm = vcl_sqrt(frobeniusNorm);
  objectnessMeasure *= 1.0 - vcl_exp(- 0.5 * vnl_math_sqr(frobeniusNorm) / vnl_math_sqr(m_Gamma));

  if (m_ScaleObjectnessMeasure) // in case, scale by the largest absolute eigenvalue
    {
    objectnessMeasure *= sortedAbsEigenValues[ImageDimension-1];
    }
  if(objectnessMeasure >= 0.0)
    {
    FinalArray.Fill(0.0);
    objit.Set( static_cast< OutputPixelType >(objectnessMeasure) ); //Here the objectness is set

    for (unsigned int j=0; j<ImageDimension; j++)
      {
      FinalArray[j] = sortedEigenVectors[j][0];
      }

    vecit.Set( static_cast< OutputArrayPixelType >(FinalArray) ); //Here are written the eigen vectors.

    RArray[0] = rA;
    RArray[1] = rB;
    RArray[2] = frobeniusNorm;
    Rit.Set(static_cast< RArrayType >(RArray));
    }
  else
    {
    FinalArray.Fill(0.0);
    objit.Set(0.0);
    vecit.Set(FinalArray);
    RArray.Fill(0.0);
    Rit.Set(static_cast< RArrayType >(RArray));
    }
  FinalArray.Fill(0);

  ++it;
  ++im;
  ++objit;
  ++vecit;
  ++Rit;
  }
}

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension >
void
EigenToObjectnessFilter< TInputPixel, TOutputPixel, VDimension >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta: " << m_Beta << std::endl;
  os << indent << "Gamma: " << m_Gamma << std::endl;
  os << indent << "ScaleObjectnessMeasure: " << m_ScaleObjectnessMeasure << std::endl;
  os << indent << "ObjectDimension: " << m_ObjectDimension << std::endl;
  os << indent << "BrightObject: " << m_BrightObject << std::endl;
}

} // end namespace itk
#endif
