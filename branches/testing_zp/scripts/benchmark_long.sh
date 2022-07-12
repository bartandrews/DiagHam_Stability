#!/bin/sh
echo "bosons 14 12 2 mono"
time ./src/Programs/QHEBosonsDeltaGround -n 2 -p 14 -l 12 > benchmark_bosons_14_12_2_mono_1.log
echo "bosons 14 12 2 bipro"
time ./src/Programs/QHEBosonsDeltaGround -S -n 2 -p 14 -l 12 > benchmark_bosons_14_12_2_bipro_1.log
echo "bosons 14 12 4 mono"
time ./src/Programs/QHEBosonsDeltaGround -n 4 -p 14 -l 12 > benchmark_bosons_14_12_4_mono_1.log
echo "bosons 14 12 4 bipro"
time ./src/Programs/QHEBosonsDeltaGround -S -n 4 -p 14 -l 12 > benchmark_bosons_14_12_4_bipro_1.log
echo "bosons 14 12 6 mono"
time ./src/Programs/QHEBosonsDeltaGround -n 6 -p 14 -l 12 > benchmark_bosons_14_12_6_mono_1.log
echo "bosons 14 12 6 bipro"
time ./src/Programs/QHEBosonsDeltaGround -S -n 6 -p 14 -l 12 > benchmark_bosons_14_12_6_bipro_1.log
echo "fermions 9 24 5 mono"
time ./src/Programs/QHEFermionsCoulomb -n 5 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_5_mono_1.log
echo "fermions 9 24 5 bipro"
time ./src/Programs/QHEFermionsCoulomb -S -n 5 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_5_bipro_1.log
echo "fermions 9 24 10 mono"
time ./src/Programs/QHEFermionsCoulomb -n 10 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_10_mono_1.log
echo "fermions 9 24 10 bipro"
time ./src/Programs/QHEFermionsCoulomb -S -n 10 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_10_bipro_1.log
echo "fermions 9 24 15 mono"
time ./src/Programs/QHEFermionsCoulomb -n 15 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_15_mono_1.log
echo "fermions 9 24 15 bipro"
time ./src/Programs/QHEFermionsCoulomb -S -n 15 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_15_bipro_1.log
echo "fermions 9 24 20 mono"
time ./src/Programs/QHEFermionsCoulomb -n 20 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_20_mono_1.log
echo "fermions 9 24 20 bipro"
time ./src/Programs/QHEFermionsCoulomb -S -n 20 -p 9 -l 24 --nbr-lz 1 > benchmark_fermions_9_24_20_bipro_1.log
