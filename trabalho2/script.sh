#!/bin/bash

all=(10 32 50 64 100 128 200 256 300 400 512 1000)
# funcoes=("Preenche", "Triangulariza")
#all=(64 100 128)
rm -rf ./TEMPO.txt ./L3.txt ./L2CACHE.txt ./FLOPS_DP.txt ./FLOPS_AVX.txt ./analise_desempenho.txt


# TEMPO
printf "\n----------- TEMPO -----------\n" >> TEMPO.txt
for i in ${all[@]} ; do
  STR1="casosTeste/entrada_${i}.in"
  echo "$STR1"
  ./calcPolinomial_tempo < $STR1 >> TEMPO.txt
  printf "\n" >> TEMPO.txt
done

echo "Terminei Tempo"

#L3
printf "\n----------- L3 -----------\n" >> L3.txt
for i in ${all[@]} ; do
  STR1="casosTeste/entrada_${i}.in"
  echo "$STR1"
  likwid-perfctr -C 3 -g L3 -m ./calcPolinomial < $STR1 | grep "L3 bandwidth" | awk -v n=$i 'BEGIN { i = 0;} { i++; print n, " L3 ", i , " - " , $6 }' >> L3.txt
  printf "\n" >> L3.txt
done

echo "Terminei L3"

#L2CACHE
printf "\n----------- L2CACHE -----------\n" >> L2CACHE.txt
for i in ${all[@]} ; do
  STR1="casosTeste/entrada_${i}.in"
  echo "$STR1"
  likwid-perfctr -C 3 -g L2CACHE -m ./calcPolinomial < $STR1 | grep "L2 miss ratio" | awk -v n=$i 'BEGIN { i = 0;} { i++; print n, " L2CACHE ", i , " - " , $6 }' >> L2CACHE.txt
  printf "\n" >> L2CACHE.txt
done

# echo "Terminei L2CACHE"

# #FLOPS_DP
printf "\n----------- FLOPS_DP -----------\n" >> FLOPS_DP.txt 
for i in ${all[@]} ; do
  STR1="casosTeste/entrada_${i}.in"
  echo "$STR1"
  likwid-perfctr -C 3 -g FLOPS_DP -m ./calcPolinomial < $STR1 | grep "|      DP MFLOP/s" | awk -v n=$i 'BEGIN { i = 0;} { i++; print n, " FLOPS_DP ", i , " - " , $5 }' >> FLOPS_DP.txt
  printf "\n" >> FLOPS_DP.txt
done

# echo "Terminei FLOPS_DP"


# #FLOPS_AVX
printf "\n----------- FLOPS_AVX -----------\n" >> FLOPS_AVX.txt 
for i in ${all[@]} ; do
  STR1="casosTeste/entrada_${i}.in"
  echo "$STR1"
  likwid-perfctr -C 3 -g FLOPS_DP -m ./calcPolinomial < $STR1 | grep "AVX DP MFLOP/s" | awk -v n=$i 'BEGIN { i = 0;} { i++; print n, " FLOPS_AVX ", i , " - " , $6 }' >> FLOPS_AVX.txt
  printf "\n" >> FLOPS_AVX.txt
done

echo "Terminei FLOPS_AVX"

cat ./TEMPO.txt ./L3.txt ./L2CACHE.txt ./FLOPS_DP.txt ./FLOPS_AVX.txt > analise_desempenho.txt

rm -rf ./TEMPO.txt ./L3.txt ./L2CACHE.txt ./FLOPS_DP.txt ./FLOPS_AVX.txt
