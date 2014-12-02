set(AngioTkFilteringTests
MultiScaleEigenObjectnessMeasureFilter.cxx
)

CreateTestDriver(AngioTkFiltering "${ITK_LIBRARIES}" "${AngioTkFilteringTests}")

angiotk_add_test(
  NAME multiScaleEigenObjectnessMeasureFilter2DTest 
  COMMAND AngioTkFilteringTestDriver
  --compare DATA{${AngiotTk_DATA_ROOT}/Baseline/ncrop2_1-1_1.mha}
            ncrop2_1-1_1.mha 
  MultiScaleEigenObjectnessMeasureFilter
  DATA{${AngiotTk_DATA_ROOT}/Input/ncrop2.tif}
  2 # dimension
  1 # sigma min
  1 # sigma max
  1 # sigma steps
  5 # gamma
  0.3 # thresholdmin
  0.4 # thresholdmax
  0 # light object
  )

angiotk_add_test(
  NAME multiScaleEigenObjectnessMeasureFilter3DTest 
  COMMAND AngioTkFilteringTestDriver
  --compare DATA{${AngiotTk_DATA_ROOT}/Baseline/cyl_1-3_4.mha}
           cyl_1-3_4.mha
  MultiScaleEigenObjectnessMeasureFilter
  DATA{${AngiotTk_DATA_ROOT}/Input/cyl.mha}
  3 # dimension
  1 # sigma min
  3 # sigma max
  4 # sigma steps
  5 # gamma
  0.2 # thresholdmin
  0.5 # thresholdmax
  25 # light object
  )
