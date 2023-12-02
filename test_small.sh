#!/bin/bash

dir="."
log="test-run.log"

N=100
m=6
eps="1e-14"

# maximal possible residual
max_re="1e-10"

if [[ $# -eq 1 ]];
then
  v="$1"
else
  v=""
fi

args1=(
  "4 4 $eps 0 ${dir}/t.ttt"
  "4 4 $eps 0 ${dir}/s.txt"
)

args2=(
  "4 4 $eps 0 ${dir}/a.txt"
  "4 4 $eps 0 ${dir}/a20.txt"
  "4 4 $eps 0 ${dir}/a40.txt"
  "4 4 $eps 0 ${dir}/b.txt"
  "6 6 $eps 0 ${dir}/c.txt"
  "6 6 $eps 0 ${dir}/d.txt"
  "6 6 $eps 0 ${dir}/e.txt"
)

max_r=$(echo "${max_re}" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

wresult_count=0
wresult_test=""
wanswer_count=0
wanswer_test=""
wapply_count=0
max_it1=0
max_it_test=""

start_time=$(date +%s)

#clear log
echo "########################################################################" | tee $log

# test on wrong files
for arg_i in "${args1[@]}"
do
  fname=$(basename $(echo "${arg_i}" | cut -d ' ' -f 5))
  echo "============================== a.out ${fname} =========================" | tee -a $log
  echo "### Run ${v} ./a.out ${arg_i}" | tee -a $log
  output=$(${v} ./a.out ${arg_i} 2>&1)
  echo "${output}" | tee -a $log
done

# test on files
for arg_i in "${args2[@]}"
do
  fname=$(basename $(echo "${arg_i}" | cut -d ' ' -f 5))
  echo "============================== a.out ${fname} =========================" | tee -a $log
  echo "### Run ${v} ./a.out ${arg_i}" | tee -a $log
  output=$(${v} ./a.out ${arg_i} 2>&1)
  echo "${output}" | tee -a $log

  result=$(echo "${output}" | grep 'a.out :')
  if [[ "${result}" == "" ]];
  then
    echo "WARNING: NO RESULT"
    (( wresult_count++ ))

    if [[ ${wresult_test} == "" ]];
    then
      wresult_test="${v} ./a.out ${arg_i}"
    fi
  else
    r1=$(echo "${result}" | cut -d ' ' -f 5 | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
    r2=$(echo "${result}" | cut -d ' ' -f 8 | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

    if [[ $(echo "${r1} == -1.0" | bc -l) -eq 1 && $(echo "${r2} == -1.0" | bc -l) -eq 1 ]];
    then
      echo "WARNING: NOT APPLICABLE"
      (( wapply_count++ ))
    else
      if [[ $(echo "${r1} < 0.0" | bc -l) -eq 1 || $(echo "${r1} >= ${max_r}" | bc -l) -eq 1
            || $(echo "${r2} < 0.0" | bc -l) -eq 1 || $(echo "${r2} >= ${max_r}" | bc -l) -eq 1 ]];
      then
        echo "WARNING: WRONG ANSWER"
        (( wanswer_count++ ))

        if [[ ${wanswer_test} == "" ]];
        then
          wanswer_test="${v} ./a.out ${arg_i}"
        fi
      fi
    fi

    it1=$(echo "${result}" | cut -d ' ' -f 14)
    if [[ $it1 -gt $max_it1 ]];
    then
      max_it1=$it1
      max_it_test="${v} ./a.out ${arg_i}"
    fi
  fi
done

# test by formula
for (( s = 1; s <= 4; s++ ))
do
  for (( n = 1; n <= $N; n++ ))
  do
    echo "============================== a.out n=$n s=$s =========================" | tee -a $log
    echo "### Run ${v} ./a.out $n $m $eps $s" | tee -a $log
    output=$(${v} ./a.out $n $m $eps $s 2>&1)
    echo "${output}" >> $log

    result=$(echo "${output}" | grep 'a.out :')
    if [[ "${result}" == "" ]];
    then
      echo "WARNING: NO RESULT"
      (( wresult_count++ ))

      if [[ ${wresult_test} == "" ]];
      then
        wresult_test="${v} ./a.out $n $m $eps $s"
      fi
    else
      echo "${result}"

      r1=$(echo "${result}" | cut -d ' ' -f 5 | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
      r2=$(echo "${result}" | cut -d ' ' -f 8 | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

      if [[ $(echo "${r1} < 0.0" | bc -l) -eq 1 || $(echo "${r1} >= ${max_r}" | bc -l) -eq 1
            || $(echo "${r2} < 0.0" | bc -l) -eq 1 || $(echo "${r2} >= ${max_r}" | bc -l) -eq 1 ]];
      then
        echo "WARNING: WRONG ANSWER"
        (( wanswer_count++ ))

        if [[ ${wanswer_test} == "" ]];
        then
          wanswer_test="${v} ./a.out $n $m $eps $s"
        fi
      fi

      it1=$(echo "${result}" | cut -d ' ' -f 14)
      if [[ $it1 -gt $max_it1 ]];
      then
        max_it1=$it1
        max_it_test="${v} ./a.out $n $m $eps $s"
      fi
    fi
  done
done

echo "########################################################################" | tee -a $log

echo "MAXIMUM ITERATIONS PER EIGENVALUE: ${max_it1} (${max_it_test})"

if [[ ${wresult_count} -gt 0 ]];
then
  echo "WARNING: POSSIBLE WRONG RESULTS: ${wresult_count} (${wresult_test})"
fi

if [[ ${wanswer_count} -gt 0 ]];
then
  echo "WARNING: POSSIBLE WRONG ANSWERS: ${wanswer_count} (${wanswer_test})"
fi

if [[ ${wapply_count} -gt 0 ]];
then
  echo "WARNING: NOT APPLICABLE: ${wapply_count} (check tests on files above)"
fi

end_time=$(date +%s)
echo "SCRIPT ELAPSED: $((end_time - start_time)) SEC"

