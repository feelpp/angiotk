/*=========================================================================

=========================================================================*/
#ifndef __itkMultiScaleEigenObjectnessFilter3D_hxx
#define __itkMultiScaleEigenObjectnessFilter3D_hxx

#include "itkMultiScaleEigenObjectnessFilter3D.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"

#define EPSILON 1e-03

namespace itk
{

/** Constructor */
template <typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
  MultiScaleEigenObjectnessFilter3D
  <TInputImage,THessianToMeasureFilter,TOutputImage>
  ::MultiScaleEigenObjectnessFilter3D()
{
  m_SigmaMin = 0.2;
  m_SigmaMax = 2.0;

  m_NumberOfSigmaSteps = 3;
  m_SigmaStepMethod = Self::LogarithmicSigmaSteps;
    
  m_Threshold = 100;
  m_UpperThreshold = NumericTraits<OutputPixelType>::Zero;
  m_LowerThreshold = NumericTraits<OutputPixelType>::Zero;
    
  // Set up neighborhood slices for all the dimensions.
  typename Neighborhood<OutputPixelType, ImageDimension>::RadiusType r;
  r.Fill(1);

  // Dummy neighborhood used to set up the slices.
  Neighborhood<OutputPixelType, ImageDimension> it;
  it.SetRadius(r);
  
  // Slice the neighborhood
  m_Center =  it.Size() / 2;
    
  //Initialize the list
  m_NodeStore = ListNodeStorageType::New();
  m_NodeList = ListType::New();

  m_HessianFilter    = HessianFilterType::New();
  m_HessianToMeasureFilter  = HessianToMeasureFilterType::New();
 
  //Instantiate Update buffer
  m_GlobalUpdateBuffer    = GlobalUpdateBufferType::New();

  this->ProcessObject::SetNumberOfRequiredOutputs(ImageDimension + 6); //Case of 3D : 9 = Ra + Rb + Rc + Vess + Magnitude + Scales + Vectors
  typename OutputImageType::Pointer vesselnessImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(0, vesselnessImage.GetPointer());

  typename OutputImageType::Pointer RAImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(1, RAImage.GetPointer());
  typename OutputImageType::Pointer RBImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(2, RBImage.GetPointer());
  typename OutputImageType::Pointer RCImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(3, RCImage.GetPointer());

  typename OutputImageType::Pointer scalesImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(4, scalesImage.GetPointer());
     
  typename OutputImageType::Pointer thresholdImage = OutputImageType::New();
  this->ProcessObject::SetNthOutput(5, thresholdImage.GetPointer());
  typename OutputImageType::Pointer eigenVector1Image = OutputImageType::New();
  this->ProcessObject::SetNthOutput(6, eigenVector1Image.GetPointer());
  typename OutputImageType::Pointer eigenVector2Image = OutputImageType::New();
  this->ProcessObject::SetNthOutput(7, eigenVector2Image.GetPointer());
  typename OutputImageType::Pointer eigenVector3Image = OutputImageType::New();
  this->ProcessObject::SetNthOutput(8, eigenVector3Image.GetPointer());

}

template <typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void
MultiScaleEigenObjectnessFilter3D
<TInputImage,THessianToMeasureFilter,TOutputImage>
::AllocateUpdateBuffer()
{
  /* The update buffer looks just like the output and holds the best response
  in the  objectness measure */

  typename TOutputImage::Pointer output = this->GetOutput();

  m_GlobalUpdateBuffer->SetSpacing(output->GetSpacing());
  m_GlobalUpdateBuffer->SetOrigin(output->GetOrigin());
  m_GlobalUpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_GlobalUpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_GlobalUpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_GlobalUpdateBuffer->Allocate();

  // Update buffer is used for > comparisons so make it really really small, just to be sure. Thanks to Hauke Heibel. 
  GlobalArrayType fillArray;
  fillArray.Fill(NumericTraits<double>::NonpositiveMin());
  m_GlobalUpdateBuffer->FillBuffer(fillArray);  
}

template <typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void
MultiScaleEigenObjectnessFilter3D
<TInputImage,THessianToMeasureFilter,TOutputImage>
::GenerateData()
{
  for(int i = 0; i < ImageDimension + 6; i++) // Allocate the output
    {
    this->GetOutput(i)->SetBufferedRegion( this->GetOutput(i)->GetRequestedRegion() );
    this->GetOutput(i)->Allocate();
    }

  AllocateUpdateBuffer(); // Allocate the buffer

  typename InputImageType::ConstPointer input = this->GetInput();

  this->m_HessianFilter->SetInput(input);
  this->m_HessianFilter->SetNormalizeAcrossScale(true);
 
  double sigma = m_SigmaMin;

  int scaleLevel = 1;

  while (sigma <= m_SigmaMax)
    {
    std::cout << "Computing measure for scale with sigma= " 
              << sigma << std::endl;
    

    m_HessianFilter->SetSigma( sigma );

    m_HessianToMeasureFilter->SetInput ( m_HessianFilter->GetOutput() ); 
    m_HessianToMeasureFilter->Update();

    this->HysteresisThresholding();
    this->UpdateMaximumResponse(sigma);
     
    sigma  = this->ComputeSigmaValue( scaleLevel );
    scaleLevel++;
     
    if ( m_NumberOfSigmaSteps == 1 )
      {
      break;
      }
    } 

  //Write out the best response to the output image
  ImageRegionIterator<GlobalUpdateBufferType> it(m_GlobalUpdateBuffer, m_GlobalUpdateBuffer->GetLargestPossibleRegion());
  it.GoToBegin();

  //Vesseleness
  ImageRegionIterator<TOutputImage> oit0(this->GetOutput(0), this->GetOutput(0)->GetLargestPossibleRegion());
  //Ra
  ImageRegionIterator<TOutputImage> oit1(this->GetOutput(1), this->GetOutput(1)->GetLargestPossibleRegion());
  //Rb
  ImageRegionIterator<TOutputImage> oit2(this->GetOutput(2), this->GetOutput(2)->GetLargestPossibleRegion());
  //Rc
  ImageRegionIterator<TOutputImage> oit3(this->GetOutput(3), this->GetOutput(3)->GetLargestPossibleRegion());
  //Scales
  //ImageRegionIterator<TOutputImage> oit4(this->GetOutput(4), this->GetOutput(4)->GetLargestPossibleRegion()); //skip scales
  //Threshold
  ImageRegionIterator<TOutputImage> oit5(this->GetOutput(5), this->GetOutput(5)->GetLargestPossibleRegion());
  //Vector1
  ImageRegionIterator<TOutputImage> oit6(this->GetOutput(6), this->GetOutput(6)->GetLargestPossibleRegion());
  //Vector2
  ImageRegionIterator<TOutputImage> oit7(this->GetOutput(7), this->GetOutput(7)->GetLargestPossibleRegion());  
  //Vector3
  ImageRegionIterator<TOutputImage> oit8(this->GetOutput(8), this->GetOutput(8)->GetLargestPossibleRegion());  
    
  oit0.GoToBegin();
  oit1.GoToBegin();
  oit2.GoToBegin();
  oit3.GoToBegin();

  oit5.GoToBegin();
  oit6.GoToBegin();
  oit7.GoToBegin();
  oit8.GoToBegin();

  GlobalArrayType globalArray;
  float objArray;
  OutputArrayType valArray, vecArray;

  while(!it.IsAtEnd())
    {
    globalArray = it.Get();

    objArray = globalArray[0];
    oit0.Set(objArray);

    oit1.Set(globalArray[1]);
    oit2.Set(globalArray[2]);
    oit3.Set(globalArray[3]);

    oit5.Set(globalArray[5]);

    for(int i = ImageDimension+5, j = ImageDimension-1; j >= 0; j--, i--)
      {
      //Get vectors starting from the last (8th output)
      vecArray[j] = globalArray[i];
      }
    oit6.Set(vecArray[0]);  
    oit7.Set(vecArray[1]);  
    oit8.Set(vecArray[2]);  
    ++oit8;
     
    ++oit0;
    ++oit1;
    ++oit2;
    ++oit3;
    
    ++oit5;
    ++oit6;
    ++oit7;
    
    ++it;
    }

  }

template <typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void
MultiScaleEigenObjectnessFilter3D
<TInputImage,THessianToMeasureFilter,TOutputImage>
::UpdateMaximumResponse(double sigma)
{
  //Global array buffer
  ImageRegionIterator<GlobalUpdateBufferType> bit(m_GlobalUpdateBuffer,m_GlobalUpdateBuffer->GetLargestPossibleRegion());
  
  //Scales output
  ImageRegionIterator<TOutputImage> oit4(this->GetOutput(4), this->GetOutput(4)->GetLargestPossibleRegion());  
  
  //Objectness measure
  ImageRegionIterator<OutputImageType> it(m_HessianToMeasureFilter->GetVesselnessOutput(),
    this->m_HessianToMeasureFilter->GetVesselnessOutput()->GetLargestPossibleRegion());

  //Vectors
  ImageRegionIterator<ArrayOutputImageType> vecit(m_HessianToMeasureFilter->GetEigenVectorsOutput(),
    this->m_HessianToMeasureFilter->GetEigenVectorsOutput()->GetLargestPossibleRegion());

  //Ra, Rb, Rc
  ImageRegionIterator<ROutputImageType> Rit(m_HessianToMeasureFilter->GetROutput(),
    this->m_HessianToMeasureFilter->GetROutput()->GetLargestPossibleRegion());
   
  //Hystheresis thresholding
  ImageRegionIterator<TOutputImage> oit5(this->GetOutput(5), this->GetOutput(5)->GetLargestPossibleRegion());

  bit.GoToBegin();
  oit4.GoToBegin();
  it.GoToBegin();
  vecit.GoToBegin();
  Rit.GoToBegin();
  oit5.GoToBegin();
   
  OutputPixelType objNew;
  OutputArrayType vectorsNew;
  RArrayType RNew;

  GlobalArrayType globalArray;

  int m, v;
  InputIndexType InputIndx;

  while(!bit.IsAtEnd())
    {
    globalArray = bit.Value();
    objNew = it.Value();
    vectorsNew = vecit.Value();  
   
    RNew = Rit.Value();
    
    InputIndx = bit.GetIndex();
    if(InputIndx[0]==48 && InputIndx[1] == 50 && InputIndx[2]==46)
      {
      std::cout << "Global array= " 
                << globalArray[0] <<" ObjNew= " << objNew << std::endl;
      }

    m = v = 0;
   
    if( globalArray[0] < objNew )
      {
      v = 1;
      globalArray[0] = objNew;

      for(int i = 1, j = 0; j < 3; i++, j++)
        {
        globalArray[i] = RNew[j];
        }
      for(int i = ImageDimension+5, j = ImageDimension-1; j >= 0; j--, i--)
        {//Get vectors starting from the last (8th output)
        globalArray[i] = vectorsNew[j]; 
        } 
      oit4.Set( sigma );
      } 
    //Hystheresis thresholding
    if(globalArray[5] < oit5.Value())
      {
      globalArray[5] = oit5.Value();
      } 
    bit.Set( globalArray );
    
    ++oit4;
    ++vecit;
    ++Rit;
    ++it;
    ++bit;  
    ++oit5;  
    }
}
     
template< typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void MultiScaleEigenObjectnessFilter3D
    <TInputImage,THessianToMeasureFilter,TOutputImage>
::HysteresisThresholding()
{
  ImageRegionIterator<TOutputImage> bit(m_HessianToMeasureFilter->GetVesselnessOutput(),
      this->m_HessianToMeasureFilter->GetVesselnessOutput()->GetLargestPossibleRegion());
  
  float value;
  ListNodeType *node;
  bit.GoToBegin();
  OutputPixelType th;
  
  while(!bit.IsAtEnd())
    {
    th = bit.Value();
    value = th;

    if(value > m_UpperThreshold)
      {
      node = m_NodeStore->Borrow();
      node->m_Value = bit.GetIndex();
      m_NodeList->PushFront(node);
      FollowEdge(bit.GetIndex());
      }
    ++bit;
    }
}

template< typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void MultiScaleEigenObjectnessFilter3D
      <TInputImage,THessianToMeasureFilter,TOutputImage>
::FollowEdge(InputIndexType index)
{
  typename OutputImageType::Pointer input = this->GetOutput(5);
   
  InputImageRegionType inputRegion = input->GetRequestedRegion();
  
  ImageRegionIteratorWithIndex<TOutputImage> bit(this->GetOutput(5), this->GetOutput(5)->GetLargestPossibleRegion());
   
  //assign iterator radius
  Size<ImageDimension> radius; 
  radius.Fill(1);
   
  ConstNeighborhoodIterator<TOutputImage> oit(radius, 
    m_HessianToMeasureFilter->GetVesselnessOutput(),
    this->m_HessianToMeasureFilter->GetVesselnessOutput()->GetLargestPossibleRegion());
  
  InputIndexType nIndex;
  InputIndexType cIndex;
  ListNodeType * node;
   
  OutputPixelType th;
  OutputPixelType nTh;

  bit.SetIndex(index);
  th = bit.Value();
   
  if(th == NumericTraits<OutputPixelType>::One )
    {
    // we must remove the node if we are not going to follow it!

    // Pop the front node from the list and read its index value.
    node = m_NodeList->Front(); // get a pointer to the first node
    m_NodeList->PopFront();  // unlink the front node
    m_NodeStore->Return(node); // return the memory for reuse
    return;
    }
  
  int nSize = m_Center * 2 + 1;  
  while(!m_NodeList->Empty())
    {
    // Pop the front node from the list and read its index value.
    node = m_NodeList->Front(); // get a pointer to the first node
    cIndex = node->m_Value;  // read the value of the first node
    m_NodeList->PopFront();  // unlink the front node
    m_NodeStore->Return(node); // return the memory for reuse

    // Move iterators to the correct index position.
    oit.SetLocation(cIndex);
    bit.SetIndex(cIndex);
    bit.Value() = 1;

    // Search the neighbors for new indicies to add to the list.
    for(int i = 0; i < nSize; i++)
      {
      nIndex = oit.GetIndex(i);
      bit.SetIndex(nIndex);
      if(inputRegion.IsInside(nIndex))
        {
        nTh = oit.GetPixel(i);
        th = bit.Value();
        if(nTh > m_LowerThreshold && th != NumericTraits<OutputPixelType>::One  )
          {
          node = m_NodeStore->Borrow();  // get a new node struct
          node->m_Value = nIndex;   // set its value
          m_NodeList->PushFront(node);   // add the new node to the list
          th = NumericTraits<OutputPixelType>::One;
          bit.SetIndex(nIndex);
          bit.Value() = th;
          }
        }
      }
    }
  }

template <typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
double
MultiScaleEigenObjectnessFilter3D
<TInputImage,THessianToMeasureFilter,TOutputImage>
::ComputeSigmaValue( int scaleLevel )
{
  double sigmaValue;
    
  if (m_NumberOfSigmaSteps < 2)
    {
    return m_SigmaMin;
    }
  
  switch (m_SigmaStepMethod)
    {
    case Self::EquispacedSigmaSteps:
      {
      const double stepSize = vnl_math_max(1e-10, (( m_SigmaMax - m_SigmaMin ) / (m_NumberOfSigmaSteps - 1)));
      sigmaValue = m_SigmaMin + stepSize * scaleLevel;
      break;
      }
    case Self::LogarithmicSigmaSteps:
      {
      const double stepSize = vnl_math_max(1e-10, (( vcl_log(m_SigmaMax) - vcl_log(m_SigmaMin) ) / (m_NumberOfSigmaSteps - 1)));
      sigmaValue = vcl_exp( vcl_log (m_SigmaMin) + stepSize * scaleLevel);
      break;
      }
    default:
      throw ExceptionObject(__FILE__, __LINE__,"Invalid SigmaStepMethod.",ITK_LOCATION);
      break;
    }
  
  return sigmaValue;
}

template<typename TInputImage, typename THessianToMeasureFilter, typename TOutputImage >
void
MultiScaleEigenObjectnessFilter3D
<TInputImage,THessianToMeasureFilter,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SigmaMin:  " << m_SigmaMin << std::endl;
  os << indent << "SigmaMax:  " << m_SigmaMax  << std::endl;
  os << indent << "NumberOfSigmaSteps:  " << m_NumberOfSigmaSteps  << std::endl;
  os << indent << "SigmaStepMethod:  " << m_SigmaStepMethod  << std::endl;
}

} // end namespace itk

#endif
