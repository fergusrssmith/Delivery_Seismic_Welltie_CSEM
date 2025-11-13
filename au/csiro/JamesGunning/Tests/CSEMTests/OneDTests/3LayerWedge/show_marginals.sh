#!/usr/bin/env bash

# plot quantile q from marginal for thickness as error bar [q,1-q], at location j, with
# true model thickness underneath. Pipe to graph -TX
#
# show_thickness_marginals.sh | graph -TX

prog='csem1DModelStudyChecker'

# 1 std dev=16
# 2std dev=2.5
quantile=2.5

j=0
for i in CMP* ; 
do
  dir="CMP"$j
  thick=`awk 'NR==7{print $1}' CMP$j/model.mod`
# thickness errors
#  $prog -Q $quantile CMP$j/ModelStudy.txt | awk -v t=$thick -v i=$j '/ m_2:/{print "#m=1,S=0\n",t,$3,"\n",t,$4; print "#m=0,S=17\n",t,t}'

# RTP errors
  $prog -M -Q $quantile CMP$j/ModelStudy.txt | awk -v t=$thick -v i=$j '/ RTP_m_2:/{print "#m=2,S=0\n",t,$3,"\n",t,$4; print "#m=0,S=17\n",t,t*100}'

  ((j++))
done



