/*=========================================================================
=========================================================================*/
#ifndef __itkEigenToObjectnessFilter_h
#define __itkEigenToObjectnessFilter_h

#include "itkFixedArray.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"

/** \class EigenToObjectnessFilter
 * \brief A filter to enhance M-dimensional objects in N-dimensional images 
 * 
 * The objectness measure is a generalization of Frangi's vesselness measure,
 * which is based on the analysis of the the Hessian eigen system. The filter
 * can enhance blob-like structures (M=0), vessel-like structures (M=1), 2D 
 * plate-like structures (M=2), hyper-plate-like structures (M=3) in N-dimensional 
 * images, with M<N.  
 * The filter takes an image of a Hessian pixels ( SymmetricSecondRankTensor pixels
 * pixels ) and produces an enhanced image. The Hessian input image can be produced 
 * using itkHessianSmoothedRecursiveGaussianImageFilter. 
 *  
 *
 * \par References
 * Frangi, AF, Niessen, WJ, Vincken, KL, & Viergever, MA (1998). Multiscale Vessel 
 * Enhancement Filtering. In Wells, WM, Colchester, A, & Delp, S, Editors, MICCAI '98 
 * Medical Image Computing and Computer-Assisted Intervention, Lecture Notes in Computer 
 * Science, pages 130-137, Springer Verlag, 1998.
 * 
 * \author Olena Tankyevych.
 * 
 * \sa HessianToObjectnessMeasureImageFilter
 * \sa MultiScaleEigenObjectnessFilter 
 * \sa MultiScaleHessianBasedMeasureImageFilter 
 * \sa Hessian3DToVesselnessMeasureImageFilter
 * \sa HessianSmoothedRecursiveGaussianImageFilter 
 * \sa SymmetricEigenAnalysisImageFilter
 * \sa SymmetricSecondRankTensor
 * 
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
namespace itk
{

template < typename TInputPixel, typename TOutputPixel, unsigned int VDimension > 
class ITK_EXPORT EigenToObjectnessFilter : public
   ImageToImageFilter< 
   Image < SymmetricSecondRankTensor< typename NumericTraits<TInputPixel>::RealType, VDimension>, VDimension>, 
   Image< FixedArray<TOutputPixel, VDimension>, VDimension > >
{
public:
  /** Standard class typedefs. */
  typedef EigenToObjectnessFilter Self;

  typedef ImageToImageFilter< Image< SymmetricSecondRankTensor
              < typename NumericTraits<TInputPixel>::RealType, VDimension>, VDimension>,
        Image< FixedArray<TOutputPixel, VDimension>, VDimension > > Superclass;

  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  typedef typename Superclass::InputImageType        InputImageType;
  typedef typename Superclass::OutputImageType       ArrayOutputImageType;
  typedef typename InputImageType::RegionType        InputImageRegionType;
  typedef typename InputImageRegionType::IndexType   InputIndexType;

  typedef typename ArrayOutputImageType::PixelType   OutputArrayPixelType;
  typedef typename InputImageType::PixelType         InputPixelType;
  typedef TOutputPixel                               OutputPixelType;  
  
  /** Image dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  typedef float                                                                     EigenValueType;
  typedef FixedArray< EigenValueType, InputPixelType::Dimension >                   EigenValueArrayType;
  typedef Image< OutputArrayPixelType, InputImageType::ImageDimension >             EigenValueImageType;
  typedef SymmetricEigenAnalysisImageFilter< InputImageType, EigenValueImageType >  EigenAnalysisFilterType;
  
  typedef FixedArray< float, 3 >                                RArrayType;
  typedef Image< RArrayType, InputImageType::ImageDimension >   ROutputImageType;
  typedef Vector< float, InputPixelType::Dimension >            EigenVectorType;
  typedef Image< EigenVectorType, ImageDimension >              EigenVectorImageType;
  typedef Image<OutputPixelType, ImageDimension>                OutputImageType;
  typedef Matrix<float, ImageDimension, ImageDimension>         AnalysisMatrixType; 
  typedef Image< AnalysisMatrixType, ImageDimension>            AnalysisMatrixImageType;
  typedef SymmetricSecondRankTensor< double, ImageDimension >   HessianTensorType;
  typedef Image< HessianTensorType >                            HessianTensorImageType;

  typedef itk::SymmetricEigenVectorAnalysisImageFilter< 
      InputImageType, 
      EigenVectorImageType,
      AnalysisMatrixImageType>   EigenVectorAnalysisFilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Get the image containing vesselness at each pixel*/
  OutputImageType* GetVesselnessOutput(void);

  /** Get the image containing the eigen vectors at each pixel*/
  ArrayOutputImageType* GetEigenVectorsOutput(void);

  /** Get the image containing RA at each pixel*/
  ROutputImageType* GetROutput(void);
  
  DataObject::Pointer MakeOutput(unsigned int idx);

  /** Set/Get Alpha, the weight corresponding to R_A (the ratio of the smallest eigenvalue that has to be large to the larger ones). 
   *  Smaller values lead to increased sensitivity to the object dimensionality. */
  itkSetMacro(Alpha,double);
  itkGetMacro(Alpha,double);

  /** Set/Get Beta, the weight corresponding to R_B (the ratio of the largest eigenvalue that has to be small to the larger ones). 
   *  Smaller values lead to increased sensitivity to the object dimensionality. */
  itkSetMacro(Beta,double);
  itkGetMacro(Beta,double);

  /** Set/Get Gamma, the weight corresponding to S (the Frobenius norm of the Hessian matrix, or second-order structureness) */
  itkSetMacro(Gamma,double);
  itkGetMacro(Gamma,double);

  /** Toggle scaling the objectness measure with the magnitude of the largest absolute eigenvalue */ 
  itkSetMacro(ScaleObjectnessMeasure,bool);
  itkGetMacro(ScaleObjectnessMeasure,bool);
  itkBooleanMacro(ScaleObjectnessMeasure);

  /** Set/Get the dimensionality of the object (0: points (blobs), 1: lines (vessels), 2: planes (plate-like structures), 3: hyper-planes. ObjectDimension must be smaller than ImageDimension. */
  itkSetMacro(ObjectDimension,int);
  itkGetMacro(ObjectDimension,int);

  /** Enhance bright structures on a dark background if true, the opposite if false. */
  itkSetMacro(BrightObject,bool);
  itkGetMacro(BrightObject,bool);

protected:
       
  EigenToObjectnessFilter();
  ~EigenToObjectnessFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Generate Data */
  void GenerateData(void);

private:

  EigenToObjectnessFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename EigenAnalysisFilterType::Pointer        m_SymmetricEigenValueFilter;
  typename EigenVectorAnalysisFilterType::Pointer  m_EigenVectorAnalysisFilter; 

  double   m_Alpha;
  double   m_Beta;
  double   m_Gamma;
  int      m_ObjectDimension;
  bool     m_BrightObject;
  bool     m_ScaleObjectnessMeasure;
  int      m_Check;
 };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEigenToObjectnessFilter.hxx"
#endif

#endif
