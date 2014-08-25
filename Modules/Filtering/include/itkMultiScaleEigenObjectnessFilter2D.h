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
#ifndef __itkMultiScaleEigenObjectnessFilter2D_h
#define __itkMultiScaleEigenObjectnessFilter2D_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkObjectStore.h"
#include "itkSparseFieldLayer.h"
#include "itkFixedArray.h"
#include <itkMinimumMaximumImageCalculator.h>

namespace itk
{

template <class TValueType>
class ListN
{
public:
  TValueType m_Value;
  ListN*     Next;
  ListN*     Previous;
};

/**\class MultiScaleEigenObjectnessFilter2D
 * \brief A filter to enhance structures using Hessian eigensystem-based 
 * measures in a multiscale framework (based on a filter "itkMultiScaleHessianBasedMeasureImageFilter")
 * 
 * The filter evaluates a Hessian-based enhancement measure, such as vesselness 
 * or objectness, at different scale levels. The Hessian-based measure is computed 
 * from the Hessian image at each scale level and the best response is selected. 
 *
 * Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
 * methods respectively. The number of scale levels is set using 
 * SetNumberOfSigmaSteps method. Exponentially distributed scale levels are 
 * computed within the bound set by the minimum and maximum sigma values 
 * 
 * The filter computes a second output image (accessed by the GetScalesOutput method)
 * containing the scales at which each pixel gave the best reponse. 
 *
 * \author Olena Tankyevych.
 *
 * \sa MultiScaleHessianBasedMeasureImageFilter
 * \sa HessianToObjectnessMeasureImageFilter 
 * \sa Hessian3DToVesselnessMeasureImageFilter 
 * \sa HessianSmoothed3DToVesselnessMeasureImageFilter 
 * \sa HessianRecursiveGaussianImageFilter 
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 * 
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
template <class TInputImage,
class THessianToMeasureFilter, 
class TOutputImage>
class ITK_EXPORT MultiScaleEigenObjectnessFilter2D 
: public
ImageToImageFilter< TInputImage, TOutputImage > 
{
public:
  
  /** Standard class typedefs. */
  typedef MultiScaleEigenObjectnessFilter2D              Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  typedef TInputImage        InputImageType;
  typedef TOutputImage       OutputImageType;

  typedef THessianToMeasureFilter           HessianToMeasureFilterType;
  typedef typename TInputImage::PixelType   InputPixelType;
  typedef typename TOutputImage::PixelType  OutputPixelType;
  
  
  typedef typename InputImageType::RegionType       InputImageRegionType;
  typedef typename InputImageRegionType::IndexType  InputIndexType;

  /** The default boundary condition is used unless overridden 
  *in the Evaluate() method. */
  typedef ZeroFluxNeumannBoundaryCondition<OutputImageType>
  DefaultBoundaryConditionType;
  typedef ConstNeighborhoodIterator<OutputImageType,
      DefaultBoundaryConditionType> NeighborhoodType;
  
  typedef ListN<InputIndexType>          ListNodeType;
  typedef ObjectStore<ListNodeType>      ListNodeStorageType;
  typedef SparseFieldLayer<ListNodeType> ListType;
  typedef typename ListType::Pointer     ListPointerType;

  /** Image dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  typedef Image<OutputPixelType, itkGetStaticConstMacro(ImageDimension)>  ScalesOutputImageType;

  /** Hessian computation filter. */
  typedef HessianRecursiveGaussianImageFilter< InputImageType >           HessianFilterType;
  typedef GradientMagnitudeRecursiveGaussianImageFilter< InputImageType >  GradientMagnitudeFilterType;
 
  typedef FixedArray<OutputPixelType, ImageDimension>       OutputArrayType; 
  typedef FixedArray<OutputPixelType, ImageDimension + 6>   GlobalArrayType; //8 = Ra + Rb + Rc + Vess + Magnitude + Scales + 2 Vectors

  typedef FixedArray< float, 3 >                               RArrayType;
  typedef Image< RArrayType, InputImageType::ImageDimension >  ROutputImageType;

  /** Update image buffer that holds the best objectness response */ 
  typedef Image< GlobalArrayType, itkGetStaticConstMacro(ImageDimension) >  GlobalUpdateBufferType;
  typedef Image< OutputArrayType, itkGetStaticConstMacro(ImageDimension) >  ArrayOutputImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get macros for SigmaMin */
  itkSetMacro(SigmaMin, double);
  itkGetMacro(SigmaMin, double);

  /** Set/Get macros for SigmaMax */
  itkSetMacro(SigmaMax, double);
  itkGetMacro(SigmaMax, double);

  /** Set/Get macros for Number of Scales */
  itkSetMacro(NumberOfSigmaSteps, int);
  itkGetMacro(NumberOfSigmaSteps, int);

  typedef enum { EquispacedSigmaSteps = 0,
    LogarithmicSigmaSteps = 1 } SigmaStepMethodType;

  /** Set/Get the method used to generate scale sequence (Equispaced or Logarithmic)*/
  itkSetMacro(SigmaStepMethod, SigmaStepMethodType);
  itkGetMacro(SigmaStepMethod, SigmaStepMethodType);

  void SetSigmaStepMethodToEquispaced()
    { 
    this->SetSigmaStepMethod(Self::EquispacedSigmaSteps);
    }
  void SetSigmaStepMethodToLogarithmic()
    { 
    this->SetSigmaStepMethod(Self::LogarithmicSigmaSteps); 
    }

  /** Get the filter used to compute the Hessian based measure */
  HessianToMeasureFilterType* GetHessianToMeasureFilter()
    {
    return m_HessianToMeasureFilter;
    }
  
  /** Get the image containing the eigen values at each pixel*/
  OutputImageType* GetVesselnessOutput()
    { 
    return  this->GetOutput(0); 
    }

  /** Set the image containing the eigen values at each pixel*/
  void SetVesselnessOutput(OutputImageType *vesselnessImage)
    { 
    this->SetNthOutput(0, vesselnessImage); 
    }

  /** Get the image containing the eigen values at each pixel*/
  OutputImageType* GetRAOutput()
   { 
   return  this->GetOutput(1); 
   }

  /** Set the image containing the eigen values at each pixel*/
  void SetRAOutput(OutputImageType *RAImage)
    { 
    this->SetNthOutput(1, RAImage); 
    }

  /** Get the image containing the eigen values at each pixel*/
  OutputImageType* GetRBOutput()
    { 
    return  this->GetOutput(2); 
    }
  /** Set the image containing the eigen values at each pixel*/
  void SetRBOutput(OutputImageType *RBImage)
    {
    this->SetNthOutput(2, RBImage); 
    }

  /** Get the image containing the eigen values at each pixel*/
  OutputImageType* GetRCOutput()
    { 
    return  this->GetOutput(3); 
    }

  /** Set the image containing the eigen values at each pixel*/
  void SetRCOutput(OutputImageType *RCImage)
    { 
    this->SetNthOutput(3, RCImage); 
    }

  /** Get the image containing the scales at which each pixel gave the best response*/
  OutputImageType* GetScalesOutput()
    { 
    return  this->GetOutput(4); 
    }

  /** Set the image containing the scales at which each pixel gave the best response*/
  void SetScalesOutput(OutputImageType *scalesImage)
    { 
    this->SetNthOutput(4, scalesImage); 
    }

  /** Get the image containing the eigen values at each pixel*/
  OutputImageType* GetMagnitudeOutput()
    { 
    return  this->GetOutput(5);
    }

  /** Set the image containing the eigen values at each pixel*/
  void SetMagnitudeOutput(OutputImageType *magnitudeImage)
    { 
    this->SetNthOutput(5, magnitudeImage); 
    }
  
  OutputImageType* GetHystheresisThresholdOutput()
    { 
    return  this->GetOutput(5); 
    }

  void SetHystheresisThresholdOutput(OutputImageType *thresholdImage)
    { 
    this->SetNthOutput(5, thresholdImage);
    }

  /** Get the image containing the eigen vectors at each pixel*/
  OutputImageType* GetEigenVector1Output()
    { 
    return  this->GetOutput(6); 
    }

  /** Set the image containing the eigen vectors at each pixel*/
  void SetEigenVector1Output(OutputImageType *eigenVector1Image)
    { 
    this->SetNthOutput(6, eigenVector1Image); 
    }

  /** Get the image containing the eigen vectors at each pixel*/
  OutputImageType* GetEigenVector2Output()
    { 
    return  this->GetOutput(7); 
    }

  /** Set the image containing the eigen vectors at each pixel*/
  void SetEigenVector2Output(OutputImageType *eigenVector2Image)
    { 
    this->SetNthOutput(7, eigenVector2Image);
    }
  
  /* Set the Threshold value for detected edges. */
  void SetThreshold(const OutputPixelType th)
    {
    this->m_Threshold = th;
    this->m_UpperThreshold = m_Threshold;
    this->m_LowerThreshold = m_Threshold/2.0;
    itkLegacyReplaceBodyMacro(SetThreshold, 2.2, SetUpperThreshold);
    }
  
  OutputPixelType GetThreshold(OutputPixelType th) 
    {
    itkLegacyReplaceBodyMacro(GetThreshold, 2.2, GetUpperThreshold);
    return this->m_Threshold; 
    }
  
  OutputPixelType GetUpperThreshold(OutputPixelType th) 
    {
    return this->m_UpperThreshold; 
    }
  
  ///* Set the Threshold value for detected edges. */
  itkSetMacro(UpperThreshold, OutputPixelType );
  itkGetConstMacro(UpperThreshold, OutputPixelType);
  
  itkSetMacro(LowerThreshold, OutputPixelType );
  itkGetConstMacro(LowerThreshold, OutputPixelType);
  
  protected:
  MultiScaleEigenObjectnessFilter2D();
  ~MultiScaleEigenObjectnessFilter2D() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Generate Data */
  void GenerateData( void );

private:

  void UpdateMaximumResponse(double sigma);
  double ComputeSigmaValue( int scaleLevel );

  void AllocateUpdateBuffer();

  /** Implement hysteresis thresholding */
  void HysteresisThresholding();
  
  /** Edge linking funciton */
  void FollowEdge(InputIndexType index);

  MultiScaleEigenObjectnessFilter2D(const Self&);   // purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double                      m_SigmaMin;
  double                      m_SigmaMax;

  int                         m_NumberOfSigmaSteps;
  SigmaStepMethodType         m_SigmaStepMethod;

  OutputPixelType m_Threshold;
  /** Upper threshold value for identifying edges. */
  OutputPixelType m_UpperThreshold;  //should be float here?
  
  /** Lower threshold value for identifying edges. */
  OutputPixelType m_LowerThreshold; //should be float here?
  
  /** "Background" value for use in thresholding. */
  OutputPixelType m_OutsideValue;

  unsigned long m_Center;
  
  typename ListNodeStorageType::Pointer           m_NodeStore;
  ListPointerType                                 m_NodeList;

  typename HessianToMeasureFilterType::Pointer    m_HessianToMeasureFilter;
  typename HessianFilterType::Pointer             m_HessianFilter;
  typename GradientMagnitudeFilterType::Pointer   m_GradientMagnitudeFilter;

  typename GlobalUpdateBufferType::Pointer        m_GlobalUpdateBuffer;
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleEigenObjectnessFilter2D.hxx"
#endif

#endif
