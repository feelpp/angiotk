TP=/home/feelpp/tp
echo "Start processing $TP..."

#SCRIPT=$TP/0_mhatonii.sh
#echo "Executing $SCRIPT"
#. $SCRIPT

#SCRIPT=$TP/1_rorpo.sh
#echo "Executing $SCRIPT"
#. $SCRIPT

#SCRIPT=$TP/1b_niitomha.sh
#echo "Executing $SCRIPT"
#. $SCRIPT

SCRIPT=$TP/2_sfi.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/3_cl.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/4_clm.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/5_ifcl.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/6_ssfi.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/7_sp.sh
echo "Executing $SCRIPT"
. $SCRIPT

SCRIPT=$TP/8_vp.sh
echo "Executing $SCRIPT"
. $SCRIPT

